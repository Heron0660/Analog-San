"""
=============================================================================
  MOTOR QUÍMICO - Analog-San
  Módulo: chem_engine.py
=============================================================================
  Contém toda a lógica de quimioinformática do projeto:

    SEÇÃO 1 - Modificações Químicas (Motor de Reações)
    SEÇÃO 2 - Cálculo de Propriedades ADMET
    SEÇÃO 3 - Filtros de Viabilidade (Lipinski, Veber)
    SEÇÃO 4 - Engine de Geração de Análogos
    SEÇÃO 5 - Exportação (CSV / SMI / MOL2) e Visualização

  -------------------------------------------------------------------------
  COMO ADICIONAR NOVAS REAÇÕES (Guia para o Pesquisador):
  -------------------------------------------------------------------------
  1. Crie uma função seguindo EXATAMENTE esta assinatura:
         def minha_reacao(mol: Chem.Mol) -> List[Chem.Mol]:
             ...
             return lista_de_mols_modificadas   # pode ser []

  2. A função recebe uma molécula RDKit e devolve uma lista de
     moléculas novas (pode ser lista vazia quando não aplicável).

  3. Adicione uma tupla ao dicionário MODIFICATION_REGISTRY:
         ("Nome Exibido", minha_reacao, (ROUTE_CLASSIC,), 3, "alta"),

  4. Pronto - o engine sorteará sua reação automaticamente.

  Opção alternativa (SMARTS):
     Use AllChem.ReactionFromSmarts(smarts_str) e chame
     rxn.RunReactants((mol,)) para obter os produtos.
     Veja os exemplos `_RXN_BENZ_TO_PYRIDINE` e
     `_RXN_BENZ_TO_PYRIMIDINE` para referência.
=============================================================================
"""

import csv
import json
import math
import random
import logging
from pathlib import Path
from typing import List, Dict, Callable, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw, rdMolDescriptors
from rdkit.Chem.rdchem import RWMol, Atom, BondType
from PIL import Image

import validation_engine as ve
from evidence_catalog import load_catalog

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
#  SEÇÃO 1: MODIFICAÇÕES QUÍMICAS (MOTOR DE REAÇÕES)
# ---------------------------------------------------------------------------

def _try_finalize(rwmol: RWMol) -> Optional[Chem.Mol]:
    """
    Tenta converter um RWMol em Mol com sanitização completa.
    Retorna None se a estrutura for quimicamente inválida.
    Encapsula exceções do RDKit de forma silenciosa.
    """
    try:
        mol = rwmol.GetMol()
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        return None


# -- 1.1  Halogenação Aromática ---------------------------------------------

def add_fluorine_aromatic(mol: Chem.Mol) -> List[Chem.Mol]:
    """
    Fluoração Aromática: insere F em cada carbono aromático com H disponível.

    Fundamento em química medicinal:
      - Aumenta a lipofilicidade (LogP +0.14 por F)
      - Bloqueia sítios de metabolismo via CYP450
      - Mimetiza ligação O-H pelo efeito polar do F
      - F é o menor substituinte que não provoca impedimento estérico
    """
    results = []
    for atom in mol.GetAtoms():
        if (atom.GetIsAromatic()
                and atom.GetAtomicNum() == 6
                and atom.GetTotalNumHs() > 0):
            rw = RWMol(mol)
            f_idx = rw.AddAtom(Atom(9))          # Flúor (Z = 9)
            rw.AddBond(atom.GetIdx(), f_idx, BondType.SINGLE)
            m = _try_finalize(rw)
            if m:
                results.append(m)
    return results


def add_chlorine_aromatic(mol: Chem.Mol) -> List[Chem.Mol]:
    """
    Cloração Aromática: insere Cl em carbonos aromáticos com H disponível.

    Fundamento em química medicinal:
      - Aumenta a lipofilicidade (LogP +0.6 por Cl)
      - Adiciona interações hidrofóbicas e halogen bonding
      - Torna o anel mais resistente à hidroxilação por CYP
    """
    results = []
    for atom in mol.GetAtoms():
        if (atom.GetIsAromatic()
                and atom.GetAtomicNum() == 6
                and atom.GetTotalNumHs() > 0):
            rw = RWMol(mol)
            cl_idx = rw.AddAtom(Atom(17))        # Cloro (Z = 17)
            rw.AddBond(atom.GetIdx(), cl_idx, BondType.SINGLE)
            m = _try_finalize(rw)
            if m:
                results.append(m)
    return results


def add_trifluoromethyl_aromatic(mol: Chem.Mol) -> List[Chem.Mol]:
    """
    Trifluorometilação Aromática: insere -CF3 em carbonos aromáticos.

    Fundamento em química medicinal:
      - Forte retirador de elétrons - modifica o pKa de grupos vizinhos
      - Aumenta substancialmente a lipofilicidade (LogP +0.9 por CF3)
      - Grupo metabolicamente estável
    """
    results = []
    for atom in mol.GetAtoms():
        if (atom.GetIsAromatic()
                and atom.GetAtomicNum() == 6
                and atom.GetTotalNumHs() > 0):
            rw = RWMol(mol)
            # Carbono central do CF3
            c_idx = rw.AddAtom(Atom(6))
            rw.AddBond(atom.GetIdx(), c_idx, BondType.SINGLE)
            # Três átomos de Flúor
            for _ in range(3):
                f_idx = rw.AddAtom(Atom(9))
                rw.AddBond(c_idx, f_idx, BondType.SINGLE)
            m = _try_finalize(rw)
            if m:
                results.append(m)
    return results


# -- 1.2  Alquilação (Metilação) --------------------------------------------

def add_methyl_aromatic(mol: Chem.Mol) -> List[Chem.Mol]:
    """
    Metilação Aromática: insere -CH3 em carbonos aromáticos com H.

    Fundamento em química medicinal:
      - Estratégia "methyl magic" - pequeno aumento de tamanho com
        grande impacto no encaixe no sítio ativo
      - Aumenta LogP moderadamente (+0.5)
      - Pode direcionar conformações bioativas
    """
    results = []
    for atom in mol.GetAtoms():
        if (atom.GetIsAromatic()
                and atom.GetAtomicNum() == 6
                and atom.GetTotalNumHs() > 0):
            rw = RWMol(mol)
            c_idx = rw.AddAtom(Atom(6))          # CH3 (carbono)
            rw.AddBond(atom.GetIdx(), c_idx, BondType.SINGLE)
            m = _try_finalize(rw)
            if m:
                results.append(m)
    return results


def replace_aromatic_halogen_with_alkyl(
    mol: Chem.Mol,
    carbon_count: int,
    branching: int = 0,
) -> List[Chem.Mol]:
    """Replace aromatic Cl, Br, or I with a small alkyl substituent.

    The transformation enumerates structures only; it does not assert that a
    particular substrate has a validated or experimentally feasible route.
    ``branching`` is 0 for a linear chain, 1 for isopropyl, and 2 for tert-butyl.
    """
    if carbon_count < 1 or branching not in (0, 1, 2):
        raise ValueError("Unsupported alkyl fragment")

    results = []
    for halogen in mol.GetAtoms():
        if halogen.GetAtomicNum() not in (17, 35, 53):
            continue
        neighbors = list(halogen.GetNeighbors())
        if len(neighbors) != 1:
            continue
        aryl_atom = neighbors[0]
        if not (aryl_atom.GetIsAromatic() and aryl_atom.GetAtomicNum() == 6):
            continue

        rw = RWMol(mol)
        halogen_idx = halogen.GetIdx()
        aryl_idx = aryl_atom.GetIdx()
        rw.RemoveAtom(halogen_idx)
        if halogen_idx < aryl_idx:
            aryl_idx -= 1

        fragment_size = 1 if branching else carbon_count
        fragment_indices = [rw.AddAtom(Atom(6)) for _ in range(fragment_size)]
        rw.AddBond(aryl_idx, fragment_indices[0], BondType.SINGLE)

        if branching == 1:
            for _ in range(2):
                methyl_idx = rw.AddAtom(Atom(6))
                rw.AddBond(fragment_indices[0], methyl_idx, BondType.SINGLE)
        elif branching == 2:
            for _ in range(3):
                methyl_idx = rw.AddAtom(Atom(6))
                rw.AddBond(fragment_indices[0], methyl_idx, BondType.SINGLE)
        else:
            for first, second in zip(fragment_indices, fragment_indices[1:]):
                rw.AddBond(first, second, BondType.SINGLE)

        candidate = _try_finalize(rw)
        if candidate:
            results.append(candidate)
    return results


def replace_aryl_halide_with_methyl(mol: Chem.Mol) -> List[Chem.Mol]:
    return replace_aromatic_halogen_with_alkyl(mol, carbon_count=1)


def replace_aryl_halide_with_ethyl(mol: Chem.Mol) -> List[Chem.Mol]:
    return replace_aromatic_halogen_with_alkyl(mol, carbon_count=2)


def replace_aryl_halide_with_isopropyl(mol: Chem.Mol) -> List[Chem.Mol]:
    return replace_aromatic_halogen_with_alkyl(mol, carbon_count=3, branching=1)


def replace_aryl_halide_with_tert_butyl(mol: Chem.Mol) -> List[Chem.Mol]:
    return replace_aromatic_halogen_with_alkyl(mol, carbon_count=4, branching=2)


def n_methylation(mol: Chem.Mol) -> List[Chem.Mol]:
    """
    N-Metilação: adiciona CH3 a átomos de nitrogênio com H disponível.

    Fundamento em química medicinal:
      - Modifica a basicidade da amina (pKa)
      - Reduz ligações de hidrogênio -> maior permeabilidade de membrana
      - Pode bloquear N-desalquilação mediada por CYP3A4
    """
    results = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetTotalNumHs() > 0:
            rw = RWMol(mol)
            c_idx = rw.AddAtom(Atom(6))
            rw.AddBond(atom.GetIdx(), c_idx, BondType.SINGLE)
            m = _try_finalize(rw)
            if m:
                results.append(m)
    return results


def o_methylation(mol: Chem.Mol) -> List[Chem.Mol]:
    """
    O-Metilação: converte -OH em -OCH3 (éter metílico).

    Fundamento em química medicinal:
      - Bloqueia glucuronidação e sulfatação de fenóis
      - Aumenta LogP (+0.5)
      - Preserva volume molecular similar ao -OH
    """
    results = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() > 0:
            rw = RWMol(mol)
            c_idx = rw.AddAtom(Atom(6))
            rw.AddBond(atom.GetIdx(), c_idx, BondType.SINGLE)
            m = _try_finalize(rw)
            if m:
                results.append(m)
    return results


# -- 1.3  Adição de Grupos Polares ------------------------------------------

def add_hydroxyl_aliphatic(mol: Chem.Mol) -> List[Chem.Mol]:
    """
    Hidroxilação Alifática: adiciona -OH a carbonos sp3 com H disponível.

    Fundamento em química medicinal:
      - Aumenta a solubilidade aquosa significativamente
      - Cria novos doadores/aceitadores de HB para interações com o alvo
      - Mimetiza produtos de metabolismo de fase I (CYP-hydroxylation)
    """
    results = []
    for atom in mol.GetAtoms():
        if (atom.GetAtomicNum() == 6
                and not atom.GetIsAromatic()
                and atom.GetTotalNumHs() > 0):
            rw = RWMol(mol)
            o_idx = rw.AddAtom(Atom(8))          # Oxigênio
            rw.AddBond(atom.GetIdx(), o_idx, BondType.SINGLE)
            m = _try_finalize(rw)
            if m:
                results.append(m)
    return results


def add_amino_aromatic(mol: Chem.Mol) -> List[Chem.Mol]:
    """
    Aminação Aromática: adiciona -NH2 a carbonos aromáticos com H.

    Fundamento em química medicinal:
      - Melhora solubilidade em pH fisiológico (protonação)
      - Adiciona doadores de HB (2 NH2 -> 2 HBDs)
      - Pode interagir com resíduos ácidos do sítio ativo (Asp, Glu)
    """
    results = []
    for atom in mol.GetAtoms():
        if (atom.GetIsAromatic()
                and atom.GetAtomicNum() == 6
                and atom.GetTotalNumHs() > 0):
            rw = RWMol(mol)
            n_idx = rw.AddAtom(Atom(7))          # Nitrogênio
            rw.AddBond(atom.GetIdx(), n_idx, BondType.SINGLE)
            m = _try_finalize(rw)
            if m:
                results.append(m)
    return results


def add_cyano_aromatic(mol: Chem.Mol) -> List[Chem.Mol]:
    """
    Cianação Aromática: adiciona -CN (nitrila, C=N) a carbonos aromáticos.

    Fundamento em química medicinal:
      - Forte retirador de elétrons; altera reatividade do anel
      - Bioisóstero de -COOH, -CONH2 e tetrazol
      - Pequeno tamanho; aumenta LogP moderadamente
    """
    results = []
    for atom in mol.GetAtoms():
        if (atom.GetIsAromatic()
                and atom.GetAtomicNum() == 6
                and atom.GetTotalNumHs() > 0):
            rw = RWMol(mol)
            c_idx = rw.AddAtom(Atom(6))          # C da nitrila
            rw.AddBond(atom.GetIdx(), c_idx, BondType.SINGLE)
            n_idx = rw.AddAtom(Atom(7))          # N = C
            rw.AddBond(c_idx, n_idx, BondType.TRIPLE)
            m = _try_finalize(rw)
            if m:
                results.append(m)
    return results


def replace_hydroxyl_with_amino(mol: Chem.Mol) -> List[Chem.Mol]:
    """
    Bioisosterismo clássico: substitui grupo hidroxila terminal por amino.

    Fundamento em química medicinal:
      - Explora a troca O/N como isosterismo clássico heteroatômico
      - Altera capacidade doadora/aceitadora de ligação de hidrogênio
      - Pode modificar basicidade, solubilidade e interação com receptor
    """
    results = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() > 0 and atom.GetDegree() == 1:
            rw = RWMol(mol)
            new_atom = rw.GetAtomWithIdx(atom.GetIdx())
            new_atom.SetAtomicNum(7)
            new_atom.SetFormalCharge(0)
            new_atom.SetNoImplicit(False)
            new_atom.SetNumExplicitHs(0)
            m = _try_finalize(rw)
            if m:
                results.append(m)
    return results


# -- 1.4  Substituições Bioisostéricas via Reações SMARTS ------------------
#
#  As reações SMARTS mapeiam átomos com (:N) para preservar a
#  conectividade com o restante da molécula durante a transformação.
#
#  Formato:  'REAGENTE >> PRODUTO'
#  Átomos mapeados garantem que as ligações externas ao fragmento
#  modificado sejam mantidas corretamente.

# Benzeno -> Piridina: substitui um -CH- aromático por N (aza-bioisóstero)
try:
    # Padrão flexível: captura qualquer cH aromático (mono/di/polissubst.)
    _RXN_BENZ_TO_PYRIDINE = AllChem.ReactionFromSmarts('[cH:1]>>[n:1]')
except Exception as e:
    _RXN_BENZ_TO_PYRIDINE = None
    logger.warning(f"Falha ao compilar reação benzeno->piridina: {e}")

# Benzeno -> Pirimidina: substitui dois -CH- por N em posições 1,3
try:
    # Substitui um cH adjacente a outro; gera padrão aza conservador
    _RXN_BENZ_TO_PYRIMIDINE = AllChem.ReactionFromSmarts('[cH:1][cH:2]>>[n:1][cH:2]')
except Exception as e:
    _RXN_BENZ_TO_PYRIMIDINE = None
    logger.warning(f"Falha ao compilar reação benzeno->pirimidina: {e}")


def _run_smarts_reaction(rxn, mol: Chem.Mol) -> List[Chem.Mol]:
    """
    Executa uma reação SMARTS sobre mol e retorna todos os produtos válidos.
    Função auxiliar usada pelas reações bioisostéricas.
    """
    if rxn is None:
        return []
    results = []
    try:
        products = rxn.RunReactants((mol,))
        for prod_tuple in products:
            for prod in prod_tuple:
                try:
                    Chem.SanitizeMol(prod)
                    results.append(prod)
                except Exception:
                    pass
    except Exception:
        pass
    return results


def benzene_to_pyridine(mol: Chem.Mol) -> List[Chem.Mol]:
    """
    Bioisosterismo: Anel benzênico -> Piridina.

    Fundamento em química medicinal:
      - Reduz LogP (mais polar) -> pode melhorar solubilidade
      - Introduz aceitador de HB adicional (N da piridina)
      - Evita metabolismo de fenóis no anel
      - Estratégia clássica de "ring nitrogen scanning"
    """
    return _run_smarts_reaction(_RXN_BENZ_TO_PYRIDINE, mol)


def benzene_to_pyrimidine(mol: Chem.Mol) -> List[Chem.Mol]:
    """
    Bioisosterismo: Anel benzênico -> Pirimidina (dois N).

    Fundamento em química medicinal:
      - Redução ainda maior do LogP vs. piridina
      - Dois aceitadores de HB - dirige orientação no sítio ativo
      - Scaffold comum em inibidores de quinases (ex: imatinibe)
    """
    return _run_smarts_reaction(_RXN_BENZ_TO_PYRIMIDINE, mol)


def ester_to_amide(mol: Chem.Mol) -> List[Chem.Mol]:
    """
    Bioisosterismo não clássico: converte éster em amida.

    Fundamento em química medicinal:
      - Aumenta estabilidade metabólica frente a esterases
      - Modifica polaridade e padrão de ligação de hidrogênio
      - Útil para modular biodisponibilidade e meia-vida
    """
    pattern = Chem.MolFromSmarts("[CX3:1](=[OX1:2])[OX2:3][#6:4]")
    results = []
    for match in mol.GetSubstructMatches(pattern):
        ester_oxygen_idx = match[2]
        rw = RWMol(mol)
        atom = rw.GetAtomWithIdx(ester_oxygen_idx)
        atom.SetAtomicNum(7)
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
        m = _try_finalize(rw)
        if m:
            results.append(m)
    return results


def ketone_to_sulfone(mol: Chem.Mol) -> List[Chem.Mol]:
    """
    Bioisosterismo não clássico: converte cetona dialquílica em sulfona.

    Fundamento em química medicinal:
      - Troca grupo carbonila por motivo sulfona mais polar e tetraédrico
      - Altera distribuição eletrônica, geometria e perfil metabólico
      - Pode melhorar seletividade por mudança conformacional
    """
    pattern = Chem.MolFromSmarts("[#6:1][CX3:2](=[OX1:3])[#6:4]")
    results = []
    seen_centers = set()
    for match in mol.GetSubstructMatches(pattern):
        carbonyl_idx = match[1]
        if carbonyl_idx in seen_centers:
            continue
        seen_centers.add(carbonyl_idx)
        rw = RWMol(mol)
        sulfur = rw.GetAtomWithIdx(carbonyl_idx)
        sulfur.SetAtomicNum(16)
        sulfur.SetFormalCharge(0)
        sulfur.SetNoImplicit(True)
        oxygen_idx = rw.AddAtom(Atom(8))
        rw.AddBond(carbonyl_idx, oxygen_idx, BondType.DOUBLE)
        m = _try_finalize(rw)
        if m:
            results.append(m)
    return results


# -- 1.5  Registro de Modificações -----------------------------------------
#
#  PARA ADICIONAR UMA NOVA REAÇÃO: basta inserir uma nova tupla abaixo.
#  Formato:
#  ("Nome", funcao, ("rota", ...), prioridade_da_regra, "nivel_manual")
#
ROUTE_CLASSIC = "classica"
ROUTE_NON_CLASSIC = "nao_classica"

MODIFICATION_REGISTRY: List[Tuple[str, Callable, Tuple[str, ...], int, str]] = [
    ("N-Metilação",                      n_methylation,                (ROUTE_NON_CLASSIC,), 1, "muito alta"),
    ("O-Metilação (-OH -> -OCH3)",       o_methylation,                (ROUTE_NON_CLASSIC,), 1, "muito alta"),
    ("Metilação Aromática (-CH3)",       add_methyl_aromatic,          (ROUTE_NON_CLASSIC,), 2, "alta"),
    ("Hidroxila->Amino (-OH -> -NH2)",   replace_hydroxyl_with_amino,  (ROUTE_CLASSIC,),     2, "alta"),
    ("Fluoração Aromática",              add_fluorine_aromatic,        (ROUTE_CLASSIC,),     3, "alta"),
    ("Hidroxilação Alifática",           add_hydroxyl_aliphatic,       (ROUTE_NON_CLASSIC,), 4, "média"),
    ("Aminação Aromática (-NH2)",        add_amino_aromatic,           (ROUTE_CLASSIC,),     4, "média"),
    ("Bioisostero: Éster->Amida",        ester_to_amide,               (ROUTE_NON_CLASSIC,), 4, "média"),
    ("Cianação Aromática (-CN)",         add_cyano_aromatic,           (ROUTE_NON_CLASSIC,), 5, "média"),
    ("Bioisostero: Benzeno->Piridina",   benzene_to_pyridine,          (ROUTE_CLASSIC,),     5, "média"),
    ("Cloração Aromática",               add_chlorine_aromatic,        (ROUTE_CLASSIC,),     5, "média"),
    ("Ar-Haleto->Metil",                  replace_aryl_halide_with_methyl,    (ROUTE_NON_CLASSIC,), 5, "média"),
    ("Ar-Haleto->Etil",                   replace_aryl_halide_with_ethyl,     (ROUTE_NON_CLASSIC,), 5, "média"),
    ("Ar-Haleto->Isopropil",              replace_aryl_halide_with_isopropyl, (ROUTE_NON_CLASSIC,), 6, "baixa"),
    ("Ar-Haleto->terc-Butil",             replace_aryl_halide_with_tert_butyl,(ROUTE_NON_CLASSIC,), 7, "baixa"),
    ("Trifluorometilação (-CF3)",        add_trifluoromethyl_aromatic, (ROUTE_CLASSIC,),     6, "baixa"),
    ("Bioisostero: Cetona->Sulfona",     ketone_to_sulfone,            (ROUTE_NON_CLASSIC,), 6, "baixa"),
    ("Bioisostero: Benzeno->Pirimidina", benzene_to_pyrimidine,        (ROUTE_CLASSIC,),     7, "baixa"),
]

_REACTION_CLASS_HYPOTHESES = {
    "N-Metilação": ("N-alquilação", "reaction_class_hypothesis"),
    "O-Metilação (-OH -> -OCH3)": ("O-alquilação / eterificação", "reaction_class_hypothesis"),
    "Metilação Aromática (-CH3)": ("Metilação C-H aromática", "reaction_class_hypothesis"),
    "Hidroxila->Amino (-OH -> -NH2)": ("Substituição O por N; rota não codificada", "graph_edit"),
    "Fluoração Aromática": ("Fluoração C-H aromática", "reaction_class_hypothesis"),
    "Hidroxilação Alifática": ("Hidroxilação C-H alifática", "reaction_class_hypothesis"),
    "Aminação Aromática (-NH2)": ("Aminação C-H aromática", "reaction_class_hypothesis"),
    "Bioisostero: Éster->Amida": ("Conversão éster-amida; rota não especificada", "reaction_class_hypothesis"),
    "Cianação Aromática (-CN)": ("Cianação C-H aromática", "reaction_class_hypothesis"),
    "Bioisostero: Benzeno->Piridina": ("Troca por aza-scaffold; rota não codificada", "scaffold_hop"),
    "Cloração Aromática": ("Cloração C-H aromática", "reaction_class_hypothesis"),
    "Ar-Haleto->Metil": ("Acoplamento aril C(sp2)-C(sp3) hipotético", "reaction_class_hypothesis"),
    "Ar-Haleto->Etil": ("Acoplamento aril C(sp2)-C(sp3) hipotético", "reaction_class_hypothesis"),
    "Ar-Haleto->Isopropil": ("Acoplamento aril C(sp2)-C(sp3) hipotético", "reaction_class_hypothesis"),
    "Ar-Haleto->terc-Butil": ("Acoplamento aril C(sp2)-C(sp3) hipotético", "reaction_class_hypothesis"),
    "Trifluorometilação (-CF3)": ("Trifluorometilação aromática", "reaction_class_hypothesis"),
    "Bioisostero: Cetona->Sulfona": ("Substituição carbonila-sulfona; rota não codificada", "graph_edit"),
    "Bioisostero: Benzeno->Pirimidina": ("Troca por diaza-scaffold; rota não codificada", "scaffold_hop"),
}

TRANSFORMATION_EVIDENCE = {
    name: {
        "status": "not_catalogued",
        "reaction_class": _REACTION_CLASS_HYPOTHESES.get(name, ("unclassified", "graph_edit"))[0],
        "record_type": _REACTION_CLASS_HYPOTHESES.get(name, ("unclassified", "graph_edit"))[1],
        "source_records": (),
        "conditions": None,
        "route_verified": False,
        "note": "Classe hipotetizada a partir da edicao estrutural; precedente experimental nao catalogado.",
    }
    for name, _, _, _, _ in MODIFICATION_REGISTRY
}


def load_transformation_evidence_catalog(path: str) -> Dict[str, Dict]:
    """Load reviewed source records without implying that they verify a route."""
    names = [name for name, _, _, _, _ in MODIFICATION_REGISTRY]
    catalog = load_catalog(path, names)
    result = {name: dict(entry) for name, entry in TRANSFORMATION_EVIDENCE.items()}
    for name, entry in catalog["transformations"].items():
        records = tuple(dict(record) for record in entry["records"])
        scopes = {record["scope"] for record in records}
        result[name].update({
            "status": "exact_transform_precedent" if "exact_transform" in scopes else (
                "reaction_class_precedent" if records else "not_catalogued"
            ),
            "source_records": records,
            "conditions": "; ".join(
                record["conditions"] for record in records if record.get("conditions")
            ) or None,
            "route_verified": False,
            "note": (
                "Precedente bibliográfico catalogado; não confirma uma rota para este substrato."
                if records else TRANSFORMATION_EVIDENCE[name]["note"]
            ),
        })
    return result


_EVIDENCE_CATALOG_PATH = Path(__file__).with_name("evidence_catalog.json")
try:
    TRANSFORMATION_EVIDENCE = load_transformation_evidence_catalog(str(_EVIDENCE_CATALOG_PATH))
except FileNotFoundError:
    pass


def get_transformation_evidence(modification_name: str) -> Dict:
    """Devolve metadados de evidência sem inferir viabilidade de síntese."""
    evidence = TRANSFORMATION_EVIDENCE.get(modification_name)
    if evidence is None:
        return {
            "status": "unknown",
            "reaction_class": None,
            "source_records": (),
            "conditions": None,
            "route_verified": False,
            "note": "Transformação sem entrada no registro de evidências.",
        }
    return dict(evidence)


def get_modification_pool(
    selected_routes: Optional[List[str]] = None,
) -> List[Tuple[str, Callable, Tuple[str, ...], int, str]]:
    """Retorna modificações habilitadas pelas rotas farmacológicas marcadas."""
    if not selected_routes:
        return sorted(MODIFICATION_REGISTRY, key=lambda item: (item[3], item[0]))

    allowed = set(selected_routes)
    return sorted([
        item for item in MODIFICATION_REGISTRY
        if allowed.intersection(item[2])
    ], key=lambda item: (item[3], item[0]))


# ---------------------------------------------------------------------------
#  SEÇÃO 2: CÁLCULO DE PROPRIEDADES ADMET
# ---------------------------------------------------------------------------

def calculate_properties(mol: Chem.Mol) -> Dict:
    """
    Calcula descritores ADMET relevantes usando o RDKit.

    Propriedades retornadas:
      MW       - Peso Molecular (Dalton)
      LogP     - Lipofilicidade (Crippen)
      HBA      - Aceitadores de ligação de hidrogênio
      HBD      - Doadores de ligação de hidrogênio
      TPSA     - Área de superfície polar topológica (A^2)
      RotBonds - Número de ligações rotativas
      HAC      - Contagem de átomos pesados (heavy atoms)
    """
    return {
        "MW":       round(Descriptors.MolWt(mol), 2),
        "LogP":     round(Descriptors.MolLogP(mol), 2),
        "HBA":      rdMolDescriptors.CalcNumHBA(mol),
        "HBD":      rdMolDescriptors.CalcNumHBD(mol),
        "TPSA":     round(rdMolDescriptors.CalcTPSA(mol), 2),
        "RotBonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "HAC":      mol.GetNumHeavyAtoms(),
    }


# ---------------------------------------------------------------------------
#  SEÇÃO 3: FILTROS DE VIABILIDADE (DRUG-LIKENESS)
# ---------------------------------------------------------------------------

def lipinski_filter(props: Dict) -> bool:
    """
    Regra dos Cinco de Lipinski (Ro5) - Lipinski et al., 1997.

    Prediz biodisponibilidade oral adequada para fármacos administrados
    por via oral (critério de violação: 0 ou 1 violação tolerada).

    Critérios obrigatórios:
      MW   <= 500 Da
      LogP <= 5
      HBD  <= 5
      HBA  <= 10
    """
    return (
        props["MW"]   <= 500 and
        props["LogP"] <= 5   and
        props["HBD"]  <= 5   and
        props["HBA"]  <= 10
    )


def veber_filter(props: Dict) -> bool:
    """
    Regras de Veber para biodisponibilidade oral - Veber et al., 2002.

    Critérios complementares a Lipinski, focados em absorção:
      RotBonds <= 10 (flexibilidade molecular)
      TPSA     <= 140 A^2 (permeabilidade de membrana)
    """
    return (
        props["RotBonds"] <= 10  and
        props["TPSA"]     <= 140
    )


def is_valid_molecule(mol) -> bool:
    """
    Verifica validade química via sanitização RDKit.

    Rejeita estruturas com:
      - Valências incorretas
      - Quiralidades impossíveis
      - Aromaticidade inconsistente
      - Qualquer outra falha de sanitização
    """
    if mol is None:
        return False
    try:
        Chem.SanitizeMol(mol)
        return True
    except Exception:
        return False


def validate_smiles(smiles: str) -> bool:
    """Verifica rapidamente se uma string SMILES é válida."""
    if not smiles or not smiles.strip():
        return False
    return Chem.MolFromSmiles(smiles.strip()) is not None


# ---------------------------------------------------------------------------
#  SEÇÃO 4: ENGINE DE GERAÇÃO DE ANÁLOGOS
# ---------------------------------------------------------------------------

def generate_analogs(
    base_smiles: str,
    target_count: int,
    use_veber: bool = True,
    selected_routes: Optional[List[str]] = None,
    progress_callback: Optional[Callable[[int, int], None]] = None,
    log_callback: Optional[Callable[[str], None]] = None,
    stop_flag: Optional[Callable[[], bool]] = None,
    random_seed: Optional[int] = ve.DEFAULT_RANDOM_SEED,
    max_sa_score: Optional[float] = 6.0,
    allow_multistep: bool = False,
    objective_profile: Optional[List[Dict[str, str]]] = None,
    selection_strategy: str = "tier_priority",
) -> List[Dict]:
    """
    Gera análogos estruturais da molécula base através de modificações químicas.

    Estratégia de geração (análogos diretos por padrão):
    ----------------------------------------------------------
      1. Inicia com a molécula base no pool de parentes.
      2. Aplica primeiro, sobre o composto-protótipo, as modificações de
         menor nível manual de prioridade da regra.
      3. Por padrão, não encadeia transformações em análogos já modificados.
      4. Aplica drug-likeness, limite opcional de SA e filtro de grupos
         explicitamente reativos antes da aceitação.
      5. Exploração multietapas é opcional e fica desligada por padrão.
      6. Controle de duplicatas via SMILES canônico (set).
      7. Ordena o resultado legado por prioridade manual da regra.

    Args:
        base_smiles:       SMILES da molécula base (será canonicalizado).
        target_count:      Número alvo de análogos aprovados.
        use_veber:         Aplica filtro de Veber além de Lipinski.
        selected_routes:   Rotas farmacológicas habilitadas na geração.
        progress_callback: fn(atual, total) - atualiza barra de progresso.
        log_callback:      fn(mensagem)     - envia mensagem ao console da GUI.
        stop_flag:         fn() -> bool      - retorna True para interromper.
        random_seed:       Seed para reprodutibilidade da exploração aleatória.
        max_sa_score:      Limite superior opcional do SA Score (escala 1-10).
        allow_multistep:   Permite combinações de modificações; padrão False.
        selection_strategy: "tier_priority" (padrão) ou "tier_diverse", que
                            seleciona por MaxMin dentro de cada tier manual.

    Returns:
        Lista de dicionários com dados de cada análogo aprovado.
        Cada dicionário tem: ID_Molecula, SMILES, Modificacao,
        PrioridadeRegra, NivelRegraManual, EvidenciaStatus, EtapaGeracao, MW, LogP,
        HBA, HBD, TPSA, RotBonds, HAC e métricas de validação científica.
    """
    # Valida e canonicaliza a molécula base
    base_mol = Chem.MolFromSmiles(base_smiles.strip())
    if base_mol is None:
        raise ValueError(f"SMILES inválido: '{base_smiles}'")

    Chem.SanitizeMol(base_mol)
    canonical_base = Chem.MolToSmiles(base_mol)
    objective_profile = objective_profile or []
    if selection_strategy not in {"tier_priority", "tier_diverse"}:
        raise ValueError("selection_strategy deve ser 'tier_priority' ou 'tier_diverse'.")
    if selection_strategy == "tier_diverse" and (allow_multistep or objective_profile):
        raise ValueError("tier_diverse requer geração direta e não pode combinar-se com frentes de Pareto.")
    if objective_profile:
        if len(objective_profile) > 2:
            raise ValueError("Defina no máximo dois objetivos moleculares.")
        ve.rank_by_objectives([], calculate_properties(base_mol), objective_profile)

    approved: List[Dict] = []
    seen_smiles = {canonical_base}
    parent_pool = [base_mol]
    max_pool_size = 300           # Evita crescimento irrestrito da memória
    modification_pool = get_modification_pool(selected_routes)
    rng = random.Random(random_seed) if random_seed is not None else random

    if not modification_pool:
        raise ValueError("Selecione pelo menos uma rota de modificação molecular.")

    # Limite de segurança para evitar loop infinito
    max_attempts = max(target_count * 120, 8000)
    attempts = 0

    def _log(msg: str):
        if log_callback:
            log_callback(msg)

    _log(f"Molécula base  : {canonical_base}")
    _log(f"Alvo           : {target_count} análogos")
    route_labels = {
        ROUTE_CLASSIC: "bioisosterismo clássico",
        ROUTE_NON_CLASSIC: "bioisosterismo não clássico",
    }
    active_routes = selected_routes or [ROUTE_CLASSIC, ROUTE_NON_CLASSIC]
    active_routes_text = ", ".join(route_labels.get(r, r) for r in active_routes)
    _log(f"Rotas ativas   : {active_routes_text}")
    _log(f"Modificações   : {len(modification_pool)} tipos disponíveis")
    if objective_profile:
        _log("Ordenação      : frentes de Pareto")
    elif selection_strategy == "tier_diverse":
        _log("Ordenação      : diversidade max-min dentro de cada tier manual (heurístico)")
    else:
        _log("Ordenação      : prioridade manual das regras (heurística)")
    _log(f"Seed           : {random_seed if random_seed is not None else 'não fixada'}")
    _log(f"Filtros        : {'Lipinski + Veber' if use_veber else 'somente Lipinski'}")
    _log(f"SA máximo      : {max_sa_score if max_sa_score is not None else 'desativado'} (heurístico; não prova síntese)")
    _log(f"Multietapas    : {'ativadas' if allow_multistep else 'desativadas; análogos diretos do protótipo'}")
    _log("-" * 58)

    def _accept_candidate(
        cand_mol: Chem.Mol,
        mod_name: str,
        rule_priority: int,
        rule_tier: str,
        generation_stage: str,
    ) -> bool:
        if not objective_profile and selection_strategy == "tier_priority" and len(approved) >= target_count:
            return False

        if not is_valid_molecule(cand_mol):
            return False

        smiles = Chem.MolToSmiles(cand_mol)
        if smiles in seen_smiles:
            return False
        seen_smiles.add(smiles)

        props = calculate_properties(cand_mol)
        passes_lipinski = lipinski_filter(props)
        passes_veber = veber_filter(props) if use_veber else True

        if not (passes_lipinski and passes_veber):
            _log(f"REJEITADO | {mod_name}: falha em Lipinski/Veber")
            return False

        sa_score = ve.synthetic_accessibility_score(cand_mol)
        if max_sa_score is not None and sa_score > max_sa_score:
            _log(f"REJEITADO | {mod_name}: SA {sa_score:.2f} > limite {max_sa_score:.2f}")
            return False

        alerts = ve.chemical_alerts(cand_mol)
        reactive_alerts = [
            alert for alert in alerts
            if any(token in alert.casefold() for token in (
                "reativo", "epóxido tensionado", "energético", "michael acceptor"
            ))
        ]
        if reactive_alerts:
            _log(f"REJEITADO | {mod_name}: alerta(s) reativo(s): {'; '.join(reactive_alerts)}")
            return False

        validation = ve.evaluate_candidate(
            base_mol=base_mol,
            analog_mol=cand_mol,
            props=props,
            rule_priority=rule_priority,
            rule_tier=rule_tier,
            use_veber=use_veber,
        )

        n = len(approved) + 1
        analog_id = f"Analog_{n:04d}"
        evidence = get_transformation_evidence(mod_name)
        analog_data = {
            "ID_Molecula": analog_id,
            "SMILES": smiles,
            "Modificacao": mod_name,
            "PrioridadeRegra": rule_priority,
            "NivelRegraManual": rule_tier,
            "EvidenciaStatus": evidence["status"],
            "EvidenciaClasseReacao": evidence["reaction_class"] or "",
            "EvidenciaTipoRegistro": evidence["record_type"],
            "EvidenciaFontesJSON": json.dumps(evidence["source_records"], ensure_ascii=False),
            "EvidenciaFonteIDs": "; ".join(
                str(record.get("source_id", "")) for record in evidence["source_records"]
                if record.get("source_id")
            ),
            "EvidenciaFonteCount": len(evidence["source_records"]),
            "EvidenciaCondicoes": evidence["conditions"] or "",
            "EvidenciaRotaVerificada": evidence["route_verified"],
            "EvidenciaObservacao": evidence["note"],
            # Deprecated aliases retained for compatibility with older exports.
            "PrioridadeSintese": rule_priority,
            "FacilidadeSintese": rule_tier,
            "EtapaGeracao": generation_stage,
            "ViabilidadeRota": "não avaliada; requer planejamento retrosintético/experimento",
            "MW": props["MW"],
            "LogP": props["LogP"],
            "HBA": props["HBA"],
            "HBD": props["HBD"],
            "TPSA": props["TPSA"],
            "RotBonds": props["RotBonds"],
            "HAC": props["HAC"],
        }
        analog_data.update(validation)
        approved.append(analog_data)

        if len(parent_pool) < max_pool_size:
            parent_pool.append(cand_mol)

        _log(
            f"OK {analog_id} | regra R{ rule_priority } {rule_tier:<10} | "
            f"Score={validation['ScoreGlobal']:>5.1f} | "
            f"{mod_name[:20]:<20} | MW={props['MW']:>6.1f}  "
            f"LogP={props['LogP']:>5.2f}"
        )

        if progress_callback and not objective_profile and selection_strategy == "tier_priority":
            progress_callback(min(len(approved), target_count), target_count)

        return True

    _log("Fase direta    : composto-protótipo por prioridade")
    for mod_name, mod_func, _, rule_priority, rule_tier in modification_pool:
        if not objective_profile and selection_strategy == "tier_priority" and len(approved) >= target_count:
            break
        if stop_flag and stop_flag():
            _log("Geração interrompida pelo usuário.")
            break

        attempts += 1
        try:
            candidates = mod_func(base_mol)
        except Exception:
            continue

        candidates = sorted(
            candidates,
            key=lambda m: Chem.MolToSmiles(m) if m is not None else "",
        )
        for cand_mol in candidates:
            if not objective_profile and selection_strategy == "tier_priority" and len(approved) >= target_count:
                break
            _accept_candidate(
                cand_mol,
                mod_name,
                rule_priority,
                rule_tier,
                "protótipo",
            )

    while allow_multistep and len(approved) < target_count and attempts < max_attempts:

        # Verifica sinal de parada enviado pela GUI
        if stop_flag and stop_flag():
            _log("Geração interrompida pelo usuário.")
            break

        attempts += 1

        # 30% das iterações partem da molécula base -> evita deriva excessiva
        if rng.random() < 0.30 or len(parent_pool) == 1:
            parent_mol = base_mol
        else:
            parent_mol = rng.choice(parent_pool)

        # Amostra regras favorecendo níveis de prioridade manual mais altos.
        mod_name, mod_func, _, rule_priority, rule_tier = rng.choices(
            modification_pool,
            weights=[1 / (item[3] ** 2) for item in modification_pool],
            k=1,
        )[0]

        # Aplica modificação -> lista de candidatos
        try:
            candidates = mod_func(parent_mol)
        except Exception:
            continue

        if not candidates:
            continue

        # Avalia cada candidato gerado pela modificação
        for cand_mol in candidates:
            if len(approved) >= target_count:
                break
            _accept_candidate(
                cand_mol,
                mod_name,
                rule_priority,
                rule_tier,
                "exploração",
            )

    _log("-" * 58)
    _log(
        f"Concluído: {len(approved)}/{target_count} análogos "
        f"em {attempts} tentativas  |  Pool final: {len(parent_pool)} mol."
    )

    if objective_profile:
        approved = ve.rank_by_objectives(approved, calculate_properties(base_mol), objective_profile)
        approved = approved[:target_count]
    elif selection_strategy == "tier_diverse":
        approved = _select_tier_diverse(approved, target_count)
        for rank, analog in enumerate(approved, start=1):
            analog["OrdemSelecaoDiversa"] = rank
        approved.sort(key=lambda analog: (
            analog["PrioridadeRegra"], analog["OrdemSelecaoDiversa"]
        ))
    else:
        approved.sort(
            key=lambda a: (
                a["PrioridadeRegra"],
                0 if a["EtapaGeracao"] == "protótipo" else 1,
                a["Modificacao"],
                a["SMILES"],
            )
        )
    for idx, analog in enumerate(approved, start=1):
        analog["ID_Molecula"] = f"Analog_{idx:04d}"

    return approved


def _select_tier_diverse(analogs: List[Dict], target_count: int) -> List[Dict]:
    """Select a diverse set within manual tiers, preserving tier precedence.

    Tiers are user-defined ordering labels, not estimates of reaction
    feasibility or synthetic difficulty.
    """
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker

    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    by_tier: Dict[int, List[Tuple[Dict, object]]] = {}
    for analog in analogs:
        mol = Chem.MolFromSmiles(analog["SMILES"])
        if mol is not None:
            by_tier.setdefault(int(analog["PrioridadeRegra"]), []).append(
                (analog, generator.GetFingerprint(mol))
            )

    selected: List[Tuple[Dict, object]] = []
    picker = MaxMinPicker()
    for tier in sorted(by_tier):
        remaining = sorted(
            by_tier[tier],
            key=lambda row: (row[0].get("Modificacao", ""), row[0].get("SMILES", "")),
        )
        slots = min(len(remaining), target_count - len(selected))
        if slots <= 0:
            break
        fingerprints = [fp for _, fp in selected] + [fp for _, fp in remaining]
        fixed_count = len(selected)
        picked_indices = picker.LazyBitVectorPick(
            fingerprints,
            len(fingerprints),
            fixed_count + slots,
            firstPicks=list(range(fixed_count)),
            seed=41242,
        )
        chosen_indices = [index - fixed_count for index in picked_indices if index >= fixed_count]
        selected_at_tier = [remaining[index] for index in chosen_indices]
        selected.extend(selected_at_tier)
        if len(selected) >= target_count:
            break
    return [analog for analog, _ in selected]


# ---------------------------------------------------------------------------
#  SEÇÃO 5: EXPORTAÇÃO E VISUALIZAÇÃO
# ---------------------------------------------------------------------------

# Colunas exportadas no CSV - compatíveis com OpenBabel / Raccoon / AutoDock
_CSV_COLUMNS = [
    "ID_Molecula", "SMILES", "Modificacao",
    "PrioridadeSintese", "FacilidadeSintese", "EtapaGeracao",
    "PrioridadeRegra", "NivelRegraManual", "EvidenciaStatus",
    "EvidenciaClasseReacao", "EvidenciaTipoRegistro", "EvidenciaFontesJSON", "EvidenciaFonteIDs",
    "EvidenciaFonteCount", "EvidenciaCondicoes", "EvidenciaRotaVerificada",
    "EvidenciaObservacao",
    "OrdemSelecaoDiversa",
    "ScoreGlobal", "CategoriaScore", "SimilaridadeTanimoto", "SAScore",
    "QtdAlertas", "AlertasQuimicos", "LipinskiOK", "VeberOK", "RegrasVioladas",
    "MW", "LogP", "HBA", "HBD", "TPSA", "RotBonds", "HAC",
    "FrentePareto", "Objetivos", "DeltaObjetivos",
]


def export_to_csv(analogs: List[Dict], filepath: str) -> None:
    """
    Exporta análogos para CSV totalmente compatível com pipelines de docking.

    Compatibilidade testada:
      - OpenBabel  : lê SMILES diretamente
      - Raccoon2   : importa CSV com SMILES + propriedades
      - AutoDock   : converte SMILES -> PDBQT via OpenBabel
      - KNIME/PyMOL: importa CSV padrão

    Uso posterior típico:
      obabel analogs.csv -O library.pdbqt --gen3d -m
    """
    if not analogs:
        raise ValueError("Nenhum análogo para exportar.")

    with open(filepath, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(analogs)


def export_to_smi(analogs: List[Dict], filepath: str) -> None:
    """
    Exporta análogos para formato .smi (SMILES canônico + ID por linha).

    Formato: <SMILES>\\t<ID>  (um composto por linha)
    Compatível com: OpenBabel, AutoDock Vina, Schrödinger Glide,
                    GOLD, DOCK, rDock e qualquer ferramenta que
                    aceite SMILES com identificadores.
    """
    if not analogs:
        raise ValueError("Nenhum análogo para exportar.")
    with open(filepath, "w", encoding="utf-8") as f:
        for a in analogs:
            f.write(f"{a['SMILES']}\t{a['ID_Molecula']}\n")


def _atom_is_amide_nitrogen(atom: Chem.Atom) -> bool:
    if atom.GetAtomicNum() != 7:
        return False
    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() != 6:
            continue
        for bond in neighbor.GetBonds():
            other = bond.GetOtherAtom(neighbor)
            if other.GetAtomicNum() == 8 and bond.GetBondType() == BondType.DOUBLE:
                return True
    return False


def _tripos_atom_type(atom: Chem.Atom) -> str:
    """Mapeamento conservador de tipos RDKit para tipos Tripos MOL2."""
    atomic_num = atom.GetAtomicNum()
    symbol = atom.GetSymbol()
    hybrid = atom.GetHybridization()

    if atomic_num == 1:
        return "H"
    if atomic_num == 6:
        if atom.GetIsAromatic():
            return "C.ar"
        if hybrid == Chem.HybridizationType.SP:
            return "C.1"
        if hybrid == Chem.HybridizationType.SP2:
            return "C.2"
        return "C.3"
    if atomic_num == 7:
        if atom.GetIsAromatic():
            return "N.ar"
        if _atom_is_amide_nitrogen(atom):
            return "N.am"
        if hybrid == Chem.HybridizationType.SP:
            return "N.1"
        if hybrid == Chem.HybridizationType.SP2:
            return "N.2"
        return "N.3"
    if atomic_num == 8:
        for bond in atom.GetBonds():
            if bond.GetBondType() == BondType.DOUBLE:
                return "O.2"
        return "O.3"
    if atomic_num == 16:
        if atom.GetIsAromatic():
            return "S.ar"
        double_oxygen = sum(
            1
            for bond in atom.GetBonds()
            if bond.GetBondType() == BondType.DOUBLE
            and bond.GetOtherAtom(atom).GetAtomicNum() == 8
        )
        if double_oxygen >= 2:
            return "S.o2"
        if double_oxygen == 1:
            return "S.o"
        return "S.3"
    if atomic_num == 15:
        return "P.3"
    if symbol in {"F", "Cl", "Br", "I"}:
        return symbol
    return symbol


def _tripos_bond_type(bond: Chem.Bond) -> str:
    if bond.GetIsAromatic():
        return "ar"
    bond_type = bond.GetBondType()
    if bond_type == BondType.SINGLE:
        return "1"
    if bond_type == BondType.DOUBLE:
        return "2"
    if bond_type == BondType.TRIPLE:
        return "3"
    return "1"


def _mol_from_smiles_with_3d(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido para MOL2: {smiles}")

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xA11A
    embed_status = AllChem.EmbedMolecule(mol, params)
    if embed_status == 0:
        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            pass
    else:
        AllChem.Compute2DCoords(mol)

    try:
        AllChem.ComputeGasteigerCharges(mol)
    except Exception:
        pass
    return mol


def _atom_charge(atom: Chem.Atom) -> float:
    if atom.HasProp("_GasteigerCharge"):
        try:
            charge = float(atom.GetProp("_GasteigerCharge"))
            if math.isfinite(charge):
                return charge
        except Exception:
            pass
    return float(atom.GetFormalCharge())


def _mol_to_mol2_block(mol: Chem.Mol, name: str) -> str:
    conf = mol.GetConformer()
    lines = [
        "@<TRIPOS>MOLECULE",
        name,
        f"{mol.GetNumAtoms():>5} {mol.GetNumBonds():>5}     1     0     0",
        "SMALL",
        "USER_CHARGES",
        "",
        "@<TRIPOS>ATOM",
    ]

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        atom_name = f"{atom.GetSymbol()}{idx + 1}"
        lines.append(
            f"{idx + 1:>7} {atom_name:<8} "
            f"{pos.x:>10.4f} {pos.y:>10.4f} {pos.z:>10.4f} "
            f"{_tripos_atom_type(atom):<6} 1 MOL {_atom_charge(atom):>10.4f}"
        )

    lines.append("@<TRIPOS>BOND")
    for bond in mol.GetBonds():
        lines.append(
            f"{bond.GetIdx() + 1:>6} "
            f"{bond.GetBeginAtomIdx() + 1:>4} "
            f"{bond.GetEndAtomIdx() + 1:>4} "
            f"{_tripos_bond_type(bond)}"
        )

    lines.extend([
        "@<TRIPOS>SUBSTRUCTURE",
        "     1 MOL         1 TEMP              0 ****  ****    0 ROOT",
        "",
    ])
    return "\n".join(lines)


def export_to_mol2(analogs: List[Dict], filepath: str) -> None:
    """
    Exporta análogos para MOL2 multi-molécula com coordenadas 3D e cargas.

    O formato é adequado para inspeção estrutural e preparação posterior em
    fluxos de docking molecular que aceitam Tripos MOL2.
    """
    if not analogs:
        raise ValueError("Nenhum análogo para exportar.")

    blocks = []
    for analog in analogs:
        smiles = analog["SMILES"]
        name = analog.get("ID_Molecula", "Analog")
        mol = _mol_from_smiles_with_3d(smiles)
        blocks.append(_mol_to_mol2_block(mol, name))

    with open(filepath, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(blocks))


def export_validation_report(
    analogs: List[Dict],
    base_smiles: str,
    output_dir: str,
    params: Optional[Dict] = None,
) -> Dict:
    """
    Gera relatório científico com ranking, CSV validado, metadados de
    reprodutibilidade e figuras PNG/SVG para artigo.
    """
    return ve.generate_scientific_report(
        analogs=analogs,
        base_smiles=base_smiles,
        output_dir=output_dir,
        params=params,
    )


def mol_to_pil_image(
    smiles: str,
    size: Tuple[int, int] = (320, 220),
    bg_color: Tuple[int, int, int] = (13, 27, 42),
) -> Optional[Image.Image]:
    """
    Renderiza um SMILES como imagem PIL 2D para exibição na GUI.

    Usa coordenadas 2D calculadas pelo RDKit (AllChem.Compute2DCoords).
    Retorna None se o SMILES for inválido ou se ocorrer qualquer erro.

    Args:
        smiles:   Código SMILES a renderizar.
        size:     Tupla (largura, altura) em pixels.
        bg_color: Cor de fundo RGB (padrão: azul escuro da GUI).
    """
    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        if mol is None:
            return None
        AllChem.Compute2DCoords(mol)
        img = Draw.MolToImage(mol, size=size, kekulize=True)
        # Aplica fundo personalizado
        bg = Image.new("RGB", size, bg_color)
        bg.paste(img, (0, 0))
        return bg
    except Exception:
        return None
