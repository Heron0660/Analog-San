"""
=============================================================================
  Analog-San - Ponto de Entrada
  Arquivo : main.py
=============================================================================
  Execute este arquivo para iniciar a aplicação:

      python main.py

  Dependências (instale via pip antes de rodar):
      pip install rdkit customtkinter pillow

  Python recomendado: 3.9 ou superior.
=============================================================================
"""

import sys


def check_dependencies() -> bool:
    """Verifica se todas as dependências obrigatórias estão instaladas."""
    missing = []
    deps = {
        "rdkit":        "rdkit",
        "customtkinter":"customtkinter",
        "PIL":          "pillow",
    }
    for module, package in deps.items():
        try:
            __import__(module)
        except ImportError:
            missing.append(package)

    if missing:
        print("=" * 60)
        print("DEPENDÊNCIAS AUSENTES")
        print("=" * 60)
        print("  Instale os pacotes abaixo e execute novamente:\n")
        print(f"  pip install {' '.join(missing)}\n")
        print("=" * 60)
        return False
    return True


def run_smoke_test(app) -> None:
    """Valida GUI, RDKit, rotas de geração e exportações CSV/MOL2."""
    import os
    import random
    import shutil
    import tempfile
    import chem_engine as ce

    mol = ce.Chem.MolFromSmiles("CCO")
    props = ce.calculate_properties(mol)
    assert props["HAC"] == 3

    img = ce.mol_to_pil_image("CCO")
    assert img is not None
    assert img.size == (320, 220)

    analog = {
        "ID_Molecula": "Smoke_0001",
        "SMILES": "CCO",
        "Modificacao": "Teste",
        "PrioridadeSintese": 1,
        "FacilidadeSintese": "muito alta",
        "EtapaGeracao": "protótipo",
        "MW": props["MW"],
        "LogP": props["LogP"],
        "HBA": props["HBA"],
        "HBD": props["HBD"],
        "TPSA": props["TPSA"],
        "RotBonds": props["RotBonds"],
        "HAC": props["HAC"],
        "ScoreGlobal": 88.5,
        "CategoriaScore": "prioritário",
        "SimilaridadeTanimoto": 1.0,
        "SAScore": 1.8,
        "QtdAlertas": 0,
        "AlertasQuimicos": "nenhum",
        "LipinskiOK": "sim",
        "VeberOK": "sim",
        "RegrasVioladas": "nenhuma",
    }
    assert ce.get_modification_pool([ce.ROUTE_CLASSIC])
    assert ce.get_modification_pool([ce.ROUTE_NON_CLASSIC])

    random.seed(1)
    generated = ce.generate_analogs(
        "c1ccccc1O",
        1,
        use_veber=True,
        selected_routes=[ce.ROUTE_CLASSIC, ce.ROUTE_NON_CLASSIC],
    )
    assert len(generated) == 1
    assert generated[0]["PrioridadeSintese"] == 1
    assert generated[0]["ScoreGlobal"] > 0
    assert "SimilaridadeTanimoto" in generated[0]
    assert "Metilação" in generated[0]["Modificacao"]

    csv_fd, csv_path = tempfile.mkstemp(suffix=".csv")
    mol2_fd, mol2_path = tempfile.mkstemp(suffix=".mol2")
    report_dir = tempfile.mkdtemp(prefix="analog_san_report_")
    os.close(csv_fd)
    os.close(mol2_fd)
    try:
        ce.export_to_csv([analog], csv_path)
        assert os.path.getsize(csv_path) > 0
        ce.export_to_mol2([analog], mol2_path)
        with open(mol2_path, "r", encoding="utf-8") as f:
            assert "@<TRIPOS>MOLECULE" in f.read()
        manifest = ce.export_validation_report(
            generated,
            "c1ccccc1O",
            report_dir,
            {"use_veber": True, "selected_routes": [ce.ROUTE_CLASSIC, ce.ROUTE_NON_CLASSIC]},
        )
        assert os.path.exists(manifest["files"]["relatorio"])
        assert os.path.exists(manifest["files"]["fig_01_png"])
    finally:
        if os.path.exists(csv_path):
            os.remove(csv_path)
        if os.path.exists(mol2_path):
            os.remove(mol2_path)
        if os.path.exists(report_dir):
            shutil.rmtree(report_dir, ignore_errors=True)

    app.withdraw()
    app.update_idletasks()
    app.destroy()


def main():
    if not check_dependencies():
        sys.exit(1)

    # Importa apenas após validação de dependências
    from gui_app import AnalogSanApp

    print("Analog-San - iniciando interface gráfica...")
    app = AnalogSanApp()
    if "--smoke-test" in sys.argv:
        run_smoke_test(app)
        print("Analog-San smoke test OK")
        return

    app.mainloop()


if __name__ == "__main__":
    main()
