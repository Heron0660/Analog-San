"""
=============================================================================
  INTERFACE GRÁFICA - Analog-San
  Módulo: gui_app.py
=============================================================================
  Interface CustomTkinter com design dark/light mode.
  Toda a lógica química está isolada em chem_engine.py.

  Arquitetura de threading:
    - Geração roda em thread daemon separada (evita travar a GUI)
    - Comunicação GUI <-> thread via queue.Queue (thread-safe)
    - GUI faz polling da queue a cada 80 ms via after()
=============================================================================
"""

import queue
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from typing import List, Dict, Optional

import customtkinter as ctk
from PIL import Image, ImageTk, ImageDraw, ImageFont

import chem_engine as ce

# -- Configuração Global do Tema --------------------------------------------
ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")


class Palette:
    """Paleta de cores centralizada - edite aqui para personalizar o tema."""
    BG_ROOT     = "#0d1117"
    BG_PANEL    = "#161b22"
    BG_CARD     = "#1c2128"
    BG_CANVAS   = "#0d1b2a"
    ACCENT_BLUE = "#1f6feb"
    ACCENT_DARK = "#0d419d"
    ACCENT_RED  = "#da3633"
    ACCENT_RED2 = "#b91c1c"
    GREEN       = "#238636"
    GREEN_HOV   = "#2ea043"
    BORDER      = "#30363d"
    TEXT        = "#e6edf3"
    TEXT_DIM    = "#7d8590"
    TEXT_GREEN  = "#56d364"
    TEXT_YELLOW = "#e3b341"
    TABLE_EVEN  = "#161b22"
    TABLE_ODD   = "#1c2128"
    SEL         = "#1f6feb"


class SafeTreeview(ttk.Treeview):
    """Treeview que evita passar -command sem valor para o Tcl/Tk."""

    def heading(self, column, option=None, **kw):
        if option == "command":
            return ""

        if kw.get("command") is None:
            kw.pop("command", None)

        return super().heading(column, option, **kw)


class AnalogSanApp(ctk.CTk):
    """
    Janela principal do Analog-San.

    Responsabilidades exclusivas desta classe:
      - Construção e layout de todos os widgets
      - Gerenciamento de estado da UI (botões habilitados/desabilitados)
      - Disparo de threads de geração
      - Polling da queue para updates thread-safe
      - Chamadas de exportação de arquivo

    Lógica química: delegada inteiramente ao módulo chem_engine.py
    """

    # ---------------------------------------------------------------------
    #  INICIALIZAÇÃO
    # ---------------------------------------------------------------------

    def __init__(self):
        super().__init__()

        self.title("Analog-San - Quimioinformática | RDKit")
        self.geometry("1340x840")
        self.minsize(1100, 720)
        self.configure(fg_color=Palette.BG_ROOT)

        # Estado interno da aplicação
        self._analogs: List[Dict] = []
        self._running: bool = False
        self._stop_requested: bool = False
        self._gui_queue: queue.Queue = queue.Queue()
        self._sort_reverse: Dict[str, bool] = {}
        self._last_base_smiles: str = ""
        self._last_generation_params: Dict = {}

        self._build_ui()
        self._poll_queue()  # Inicia loop de polling (seguro para threads)

    # ---------------------------------------------------------------------
    #  CONSTRUÇÃO DA INTERFACE
    # ---------------------------------------------------------------------

    def _build_ui(self):
        """Ponto de entrada para construção de todos os componentes."""
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)
        self.grid_rowconfigure(2, weight=0)

        self._build_header()
        self._build_body()
        self._build_statusbar()

    # -- Header ------------------------------------------------------------

    def _build_header(self):
        hdr = ctk.CTkFrame(
            self, fg_color=Palette.BG_PANEL, corner_radius=0,
            border_width=1, border_color=Palette.BORDER, height=56,
        )
        hdr.grid(row=0, column=0, sticky="ew")
        hdr.grid_propagate(False)
        hdr.grid_columnconfigure(1, weight=1)

        # Logo / título
        ctk.CTkLabel(
            hdr,
            text="Analog-San",
            font=ctk.CTkFont(family="Segoe UI", size=20, weight="bold"),
            text_color=Palette.TEXT,
        ).grid(row=0, column=0, padx=20, pady=10, sticky="w")

        ctk.CTkLabel(
            hdr,
            text="Plataforma de Quimioinformática para Otimização de Candidatos a Fármaco",
            font=ctk.CTkFont(size=11),
            text_color=Palette.TEXT_DIM,
        ).grid(row=0, column=1, padx=10, sticky="w")

        # Botão de tema
        self._theme_btn = ctk.CTkButton(
            hdr, text="Modo Claro", width=110, height=30,
            font=ctk.CTkFont(size=11),
            fg_color=Palette.BG_CARD, hover_color=Palette.BORDER,
            border_width=1, border_color=Palette.BORDER,
            command=self._toggle_theme,
        )
        self._theme_btn.grid(row=0, column=2, padx=12, pady=10)

    # -- Body (dois painéis) -----------------------------------------------

    def _build_body(self):
        body = ctk.CTkFrame(self, fg_color="transparent")
        body.grid(row=1, column=0, sticky="nsew", padx=10, pady=(8, 4))
        body.grid_columnconfigure(0, weight=0, minsize=350)
        body.grid_columnconfigure(1, weight=1)
        body.grid_rowconfigure(0, weight=1)

        self._build_left_panel(body)
        self._build_right_panel(body)

    # -- Painel Esquerdo: Controles ----------------------------------------

    def _build_left_panel(self, parent):
        frame = ctk.CTkScrollableFrame(
            parent, fg_color=Palette.BG_PANEL,
            scrollbar_button_color=Palette.BORDER,
            corner_radius=10, width=340,
            border_width=1, border_color=Palette.BORDER,
        )
        frame.grid(row=0, column=0, sticky="nsew", padx=(0, 8))
        frame.grid_columnconfigure(0, weight=1)

        row = 0

        # -- Entrada SMILES ------------------------------------------------
        row = self._section_header(frame, "Molécula Base (SMILES)", row)

        self._smiles_entry = ctk.CTkEntry(
            frame, placeholder_text="Cole o código SMILES aqui...",
            height=40, font=ctk.CTkFont(family="Consolas", size=12),
            fg_color=Palette.BG_CARD, border_color=Palette.BORDER,
            text_color=Palette.TEXT,
        )
        self._smiles_entry.grid(row=row, column=0, padx=12, pady=(4, 6), sticky="ew")
        self._smiles_entry.bind("<KeyRelease>", self._on_smiles_change)
        row += 1

        # -- Botões de moléculas exemplo -----------------------------------
        row = self._section_header(frame, "Exemplos Rápidos", row)

        eg_frame = ctk.CTkFrame(frame, fg_color="transparent")
        eg_frame.grid(row=row, column=0, sticky="ew", padx=12, pady=(4, 8))
        eg_frame.grid_columnconfigure((0, 1, 2), weight=1)
        row += 1

        EXAMPLES = [
            ("Aspirina",    "CC(=O)Oc1ccccc1C(=O)O"),
            ("Ibuprofeno",  "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
            ("Cafeína",     "Cn1cnc2c1c(=O)n(c(=O)n2C)C"),
            ("Metformina",  "CN(C)C(=N)NC(=N)N"),
            ("Paracetamol", "CC(=O)Nc1ccc(O)cc1"),
            ("Atorvastatina","CC(C)c1c(C(=O)Nc2ccccc2F)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O"),
        ]
        for i, (name, smi) in enumerate(EXAMPLES):
            col = i % 3
            ctk.CTkButton(
                eg_frame, text=name, height=28,
                font=ctk.CTkFont(size=10),
                fg_color=Palette.BG_CARD, hover_color=Palette.ACCENT_DARK,
                border_width=1, border_color=Palette.BORDER,
                command=lambda s=smi: self._load_example(s),
            ).grid(row=i // 3, column=col, padx=3, pady=3, sticky="ew")

        # -- Visualizador 2D -----------------------------------------------
        row = self._section_header(frame, "Visualização 2D", row)

        canvas_frame = ctk.CTkFrame(
            frame, fg_color=Palette.BG_CANVAS,
            border_width=1, border_color=Palette.BORDER, corner_radius=8,
        )
        canvas_frame.grid(row=row, column=0, padx=12, pady=(4, 10), sticky="ew")
        row += 1

        self._mol_canvas = tk.Canvas(
            canvas_frame, width=316, height=210,
            bg=Palette.BG_CANVAS, highlightthickness=0,
        )
        self._mol_canvas.pack(padx=2, pady=2)
        self._canvas_img_ref = None  # Referência para evitar coleta pelo GC

        self._canvas_hint()

        # -- Parâmetros de Geração -----------------------------------------
        row = self._section_header(frame, "Parâmetros de Geração", row)

        # Slider de quantidade
        qty_lbl_frame = ctk.CTkFrame(frame, fg_color="transparent")
        qty_lbl_frame.grid(row=row, column=0, padx=12, sticky="ew")
        qty_lbl_frame.grid_columnconfigure(1, weight=1)
        row += 1

        ctk.CTkLabel(
            qty_lbl_frame, text="No. de Análogos:",
            font=ctk.CTkFont(size=12), text_color=Palette.TEXT_DIM,
        ).grid(row=0, column=0, sticky="w")

        self._qty_val_lbl = ctk.CTkLabel(
            qty_lbl_frame, text="100",
            font=ctk.CTkFont(size=14, weight="bold"),
            text_color=Palette.ACCENT_BLUE,
        )
        self._qty_val_lbl.grid(row=0, column=2, padx=(8, 0))

        self._qty_slider = ctk.CTkSlider(
            frame, from_=1, to=1000, number_of_steps=999,
            command=self._on_qty_slider,
            button_color=Palette.ACCENT_BLUE,
            progress_color=Palette.ACCENT_BLUE,
        )
        self._qty_slider.set(100)
        self._qty_slider.grid(row=row, column=0, padx=12, pady=(4, 8), sticky="ew")
        row += 1

        # Entrada numérica precisa
        spin_frame = ctk.CTkFrame(frame, fg_color="transparent")
        spin_frame.grid(row=row, column=0, padx=12, pady=(0, 10), sticky="ew")
        spin_frame.grid_columnconfigure(1, weight=1)
        row += 1

        ctk.CTkLabel(
            spin_frame, text="Valor exato (1-1000):",
            font=ctk.CTkFont(size=11), text_color=Palette.TEXT_DIM,
        ).grid(row=0, column=0, sticky="w")

        self._qty_entry = ctk.CTkEntry(
            spin_frame, width=80, height=30,
            font=ctk.CTkFont(family="Consolas", size=12),
            fg_color=Palette.BG_CARD, border_color=Palette.BORDER,
        )
        self._qty_entry.insert(0, "100")
        self._qty_entry.grid(row=0, column=1, padx=(8, 0), sticky="e")
        self._qty_entry.bind("<Return>",       self._on_qty_entry_enter)
        self._qty_entry.bind("<FocusOut>",     self._on_qty_entry_enter)

        # -- Filtros -------------------------------------------------------
        row = self._section_header(frame, "Filtros de Drug-Likeness", row)

        filter_frame = ctk.CTkFrame(frame, fg_color=Palette.BG_CARD,
                                    corner_radius=8, border_width=1,
                                    border_color=Palette.BORDER)
        filter_frame.grid(row=row, column=0, padx=12, pady=(4, 10), sticky="ew")
        row += 1

        self._lipinski_var = ctk.BooleanVar(value=True)
        self._veber_var    = ctk.BooleanVar(value=True)
        self._classic_route_var = ctk.BooleanVar(value=True)
        self._non_classic_route_var = ctk.BooleanVar(value=True)
        self._multistep_var = ctk.BooleanVar(value=False)
        self._max_sa_var = ctk.StringVar(value="6.0")

        ctk.CTkCheckBox(
            filter_frame,
            text="Lipinski Ro5  (MW<=500, LogP<=5, HBD<=5, HBA<=10)",
            variable=self._lipinski_var, state="disabled",
            font=ctk.CTkFont(size=11), checkmark_color="white",
            fg_color=Palette.ACCENT_BLUE,
        ).pack(anchor="w", padx=10, pady=(8, 4))

        ctk.CTkCheckBox(
            filter_frame,
            text="Veber  (RotBonds<=10, TPSA<=140 A^2)",
            variable=self._veber_var,
            font=ctk.CTkFont(size=11), checkmark_color="white",
            fg_color=Palette.ACCENT_BLUE,
        ).pack(anchor="w", padx=10, pady=(0, 4))

        ctk.CTkCheckBox(
            filter_frame,
            text="Rota clássica (bioisósteres clássicos)",
            variable=self._classic_route_var,
            font=ctk.CTkFont(size=11), checkmark_color="white",
            fg_color=Palette.ACCENT_BLUE,
        ).pack(anchor="w", padx=10, pady=(0, 4))

        ctk.CTkCheckBox(
            filter_frame,
            text="Rota não clássica (bioisósteres funcionais)",
            variable=self._non_classic_route_var,
            font=ctk.CTkFont(size=11), checkmark_color="white",
            fg_color=Palette.ACCENT_BLUE,
        ).pack(anchor="w", padx=10, pady=(0, 8))

        quality_frame = ctk.CTkFrame(filter_frame, fg_color="transparent")
        quality_frame.pack(fill="x", padx=10, pady=(0, 8))
        ctk.CTkLabel(
            quality_frame,
            text="SA Score máximo (1-10; estimativa, não comprova síntese):",
            font=ctk.CTkFont(size=11),
        ).pack(side="left")
        ctk.CTkEntry(
            quality_frame, textvariable=self._max_sa_var, width=58, height=28,
        ).pack(side="left", padx=(8, 0))
        ctk.CTkCheckBox(
            filter_frame,
            text="Permitir combinações de modificações (exploração multietapas)",
            variable=self._multistep_var,
            font=ctk.CTkFont(size=11), checkmark_color="white",
            fg_color=Palette.ACCENT_BLUE,
        ).pack(anchor="w", padx=10, pady=(0, 8))

        row = self._section_header(frame, "Objetivos de otimização (opcional)", row)
        self._objective_vars = []
        objective_frame = ctk.CTkFrame(frame, fg_color=Palette.BG_CARD,
                                       corner_radius=8, border_width=1,
                                       border_color=Palette.BORDER)
        objective_frame.grid(row=row, column=0, padx=12, pady=(4, 10), sticky="ew")
        row += 1
        prop_values = ["Nenhum", "LogP", "MW", "TPSA", "HBD", "HBA", "RotBonds"]
        direction_values = ["Aumentar", "Reduzir", "Manter"]
        for index in range(2):
            prop_var = ctk.StringVar(value="Nenhum")
            direction_var = ctk.StringVar(value="Aumentar")
            self._objective_vars.append((prop_var, direction_var))
            ctk.CTkOptionMenu(
                objective_frame, values=prop_values, variable=prop_var,
                width=105, height=28, font=ctk.CTkFont(size=11),
            ).grid(row=index, column=0, padx=(8, 4), pady=5, sticky="ew")
            ctk.CTkOptionMenu(
                objective_frame, values=direction_values, variable=direction_var,
                width=105, height=28, font=ctk.CTkFont(size=11),
            ).grid(row=index, column=1, padx=(4, 8), pady=5, sticky="ew")
        objective_frame.grid_columnconfigure((0, 1), weight=1)

        # -- Botões de Ação ------------------------------------------------
        row = self._section_header(frame, "Ações", row)

        self._gen_btn = ctk.CTkButton(
            frame,
            text="Gerar Análogos",
            height=46, font=ctk.CTkFont(size=14, weight="bold"),
            fg_color=Palette.GREEN, hover_color=Palette.GREEN_HOV,
            command=self._on_generate,
        )
        self._gen_btn.grid(row=row, column=0, padx=12, pady=(4, 6), sticky="ew")
        row += 1

        action2_frame = ctk.CTkFrame(frame, fg_color="transparent")
        action2_frame.grid(row=row, column=0, padx=12, pady=(0, 12), sticky="ew")
        action2_frame.grid_columnconfigure((0, 1, 2), weight=1)
        row += 1

        self._stop_btn = ctk.CTkButton(
            action2_frame, text="Parar", height=36,
            font=ctk.CTkFont(size=12),
            fg_color=Palette.ACCENT_RED, hover_color=Palette.ACCENT_RED2,
            state="disabled", command=self._on_stop,
        )
        self._stop_btn.grid(row=0, column=0, sticky="ew", padx=(0, 4))

        self._export_btn = ctk.CTkButton(
            action2_frame, text="Exportar", height=36,
            font=ctk.CTkFont(size=12),
            fg_color=Palette.ACCENT_DARK, hover_color=Palette.ACCENT_BLUE,
            state="disabled", command=self._on_export,
        )
        self._export_btn.grid(row=0, column=1, sticky="ew", padx=4)

        self._report_btn = ctk.CTkButton(
            action2_frame, text="Relatório", height=36,
            font=ctk.CTkFont(size=12),
            fg_color=Palette.BG_CARD, hover_color=Palette.BORDER,
            border_width=1, border_color=Palette.BORDER,
            state="disabled", command=self._on_report,
        )
        self._report_btn.grid(row=0, column=2, sticky="ew", padx=(4, 0))

        # -- Barra de Progresso --------------------------------------------
        row = self._section_header(frame, "Progresso", row)

        self._progress_bar = ctk.CTkProgressBar(
            frame, height=16, corner_radius=8,
            progress_color=Palette.ACCENT_BLUE,
            fg_color=Palette.BG_CARD,
        )
        self._progress_bar.set(0)
        self._progress_bar.grid(row=row, column=0, padx=12, pady=(4, 4), sticky="ew")
        row += 1

        self._progress_lbl = ctk.CTkLabel(
            frame, text="Aguardando...",
            font=ctk.CTkFont(size=11), text_color=Palette.TEXT_DIM,
        )
        self._progress_lbl.grid(row=row, column=0, pady=(0, 16))
        row += 1

    # -- Painel Direito: Resultados ----------------------------------------

    def _build_right_panel(self, parent):
        frame = ctk.CTkFrame(
            parent, fg_color=Palette.BG_PANEL,
            corner_radius=10, border_width=1, border_color=Palette.BORDER,
        )
        frame.grid(row=0, column=1, sticky="nsew")
        frame.grid_columnconfigure(0, weight=1)
        frame.grid_rowconfigure(1, weight=3)
        frame.grid_rowconfigure(3, weight=1)

        # -- Cabeçalho da tabela + estatísticas ---------------------------
        tbl_hdr = ctk.CTkFrame(frame, fg_color="transparent")
        tbl_hdr.grid(row=0, column=0, sticky="ew", padx=12, pady=(12, 4))
        tbl_hdr.grid_columnconfigure(1, weight=1)

        ctk.CTkLabel(
            tbl_hdr, text="Análogos Gerados",
            font=ctk.CTkFont(size=13, weight="bold"),
            text_color=Palette.TEXT_DIM,
        ).grid(row=0, column=0, sticky="w")

        self._stats_lbl = ctk.CTkLabel(
            tbl_hdr, text="",
            font=ctk.CTkFont(size=11), text_color=Palette.TEXT_YELLOW,
        )
        self._stats_lbl.grid(row=0, column=1, sticky="e")

        # -- Treeview (Tabela de Resultados) ------------------------------
        tbl_container = ctk.CTkFrame(
            frame, fg_color=Palette.BG_CANVAS,
            corner_radius=8, border_width=1, border_color=Palette.BORDER,
        )
        tbl_container.grid(row=1, column=0, sticky="nsew", padx=12, pady=(0, 8))
        tbl_container.grid_columnconfigure(0, weight=1)
        tbl_container.grid_rowconfigure(0, weight=1)

        self._configure_treeview_style()

        columns = (
            "ID", "SMILES", "Modificacao", "Regra", "Pareto", "Objetivo",
            "Score", "Tanimoto", "SA", "Alertas",
            "MW", "LogP", "HBA", "HBD", "TPSA", "RotBonds",
        )
        self._table = SafeTreeview(
            tbl_container, columns=columns, show="headings",
            style="AnalogSan.Treeview", selectmode="browse",
        )

        # Configuração das colunas - (header_text, width, anchor, sortable)
        col_cfg = {
            "ID":         ("ID",           95,  "center", True),
            "SMILES":     ("SMILES",       290, "w",      False),
            "Modificacao":("Modificação",  155, "w",      True),
            "Regra":      ("Prior. regra",  105, "center", True),
            "Pareto":     ("Pareto",        60, "center", True),
            "Objetivo":   ("Δ objetivo",    150, "w",      False),
            "Score":      ("Score (heur.)", 88, "center", True),
            "Tanimoto":   ("Tanimoto",      82, "center", True),
            "SA":         ("SA",            58, "center", True),
            "Alertas":    ("Alertas",       70, "center", True),
            "MW":         ("MW (Da)",       82, "center", True),
            "LogP":       ("LogP",          72, "center", True),
            "HBA":        ("HBA",           56, "center", True),
            "HBD":        ("HBD",           56, "center", True),
            "TPSA":       ("TPSA (A^2)",    85, "center", True),
            "RotBonds":   ("Rot.",          52, "center", True),
        }
        self._sortable_columns = {
            col for col, (_, _, _, sortable) in col_cfg.items() if sortable
        }
        for col, (header, width, anchor, sortable) in col_cfg.items():
            heading_text = header + (" " if sortable else "")
            self._table.heading(col, text=heading_text)
            self._table.column(col, width=width, anchor=anchor, minwidth=40,
                                stretch=(col == "SMILES"))

        self._table.tag_configure("even", background="#1c2128")
        self._table.tag_configure("odd",  background="#161b22")
        self._table.tag_configure("sel",  background=Palette.SEL)

        vsb = ttk.Scrollbar(tbl_container, orient="vertical",
                             command=self._table.yview)
        hsb = ttk.Scrollbar(tbl_container, orient="horizontal",
                             command=self._table.xview)
        self._table.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        self._table.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, columnspan=2, sticky="ew")

        self._table.bind("<<TreeviewSelect>>", self._on_row_select)
        self._table.bind("<Button-1>", self._on_table_click, add="+")

        # -- Detalhe da molécula selecionada ------------------------------
        detail_hdr = ctk.CTkFrame(frame, fg_color="transparent")
        detail_hdr.grid(row=2, column=0, sticky="ew", padx=12, pady=(0, 4))

        ctk.CTkLabel(
            detail_hdr, text="Detalhe - Molécula Selecionada",
            font=ctk.CTkFont(size=13, weight="bold"),
            text_color=Palette.TEXT_DIM,
        ).pack(side="left")

        detail_frame = ctk.CTkFrame(
            frame, fg_color=Palette.BG_CARD,
            corner_radius=8, border_width=1, border_color=Palette.BORDER,
        )
        detail_frame.grid(row=3, column=0, sticky="nsew", padx=12, pady=(0, 8))
        detail_frame.grid_columnconfigure(0, weight=1)
        detail_frame.grid_columnconfigure(1, weight=0)
        detail_frame.grid_rowconfigure(0, weight=1)

        self._detail_lbl = ctk.CTkLabel(
            detail_frame,
            text="<- Selecione uma linha na tabela para ver detalhes",
            font=ctk.CTkFont(family="Consolas", size=11),
            text_color=Palette.TEXT_DIM, justify="left", anchor="nw",
        )
        self._detail_lbl.grid(row=0, column=0, padx=14, pady=10, sticky="nw")

        self._detail_canvas = tk.Canvas(
            detail_frame, width=220, height=130,
            bg=Palette.BG_CANVAS, highlightthickness=0,
        )
        self._detail_canvas.grid(row=0, column=1, padx=8, pady=6, sticky="e")
        self._detail_img_ref = None

        # -- Console de Log ------------------------------------------------
        log_hdr = ctk.CTkFrame(frame, fg_color="transparent")
        log_hdr.grid(row=4, column=0, sticky="ew", padx=12, pady=(0, 4))
        ctk.CTkLabel(
            log_hdr, text="Console de Log",
            font=ctk.CTkFont(size=13, weight="bold"),
            text_color=Palette.TEXT_DIM,
        ).pack(side="left")

        ctk.CTkButton(
            log_hdr, text="Limpar", width=80, height=24,
            font=ctk.CTkFont(size=10),
            fg_color=Palette.BG_CARD, hover_color=Palette.BORDER,
            command=self._clear_log,
        ).pack(side="right")

        self._log_box = ctk.CTkTextbox(
            frame, height=130,
            font=ctk.CTkFont(family="Consolas", size=11),
            fg_color=Palette.BG_CANVAS, text_color=Palette.TEXT_GREEN,
            border_width=1, border_color=Palette.BORDER,
            state="disabled", wrap="none",
        )
        self._log_box.grid(row=5, column=0, sticky="ew", padx=12, pady=(0, 12))

    # -- Status Bar --------------------------------------------------------

    def _build_statusbar(self):
        sb = ctk.CTkFrame(
            self, fg_color=Palette.BG_PANEL, corner_radius=0,
            border_width=1, border_color=Palette.BORDER, height=26,
        )
        sb.grid(row=2, column=0, sticky="ew")
        sb.grid_propagate(False)
        sb.grid_columnconfigure(1, weight=1)

        ctk.CTkLabel(
            sb,
            text="Analog-San  |  RDKit + CustomTkinter  "
                 "|  Ordem: prioridade manual da regra  "
                 "|  Exporta: .csv / .smi / .mol2 / relatório científico",
            font=ctk.CTkFont(size=9), text_color=Palette.TEXT_DIM,
        ).grid(row=0, column=0, padx=10, sticky="w")

        ctk.CTkLabel(
            sb,
            text="By Heron",
            font=ctk.CTkFont(size=9, weight="bold"),
            text_color=Palette.TEXT_YELLOW,
        ).grid(row=0, column=2, padx=10, sticky="e")

    # ---------------------------------------------------------------------
    #  HELPERS DE LAYOUT
    # ---------------------------------------------------------------------

    def _section_header(self, parent, text: str, row: int) -> int:
        """Insere um rótulo de seção estilizado e retorna a próxima linha."""
        ctk.CTkLabel(
            parent, text=text, anchor="w",
            font=ctk.CTkFont(size=11, weight="bold"),
            text_color=Palette.TEXT_DIM,
        ).grid(row=row, column=0, padx=12, pady=(14, 2), sticky="ew")
        return row + 1

    def _configure_treeview_style(self):
        """Estiliza o ttk.Treeview para combinar com o tema escuro."""
        style = ttk.Style()
        style.theme_use("clam")
        style.configure(
            "AnalogSan.Treeview",
            background=Palette.BG_CANVAS,
            foreground=Palette.TEXT,
            rowheight=26,
            fieldbackground=Palette.BG_CANVAS,
            borderwidth=0,
            font=("Consolas", 10),
        )
        style.configure(
            "AnalogSan.Treeview.Heading",
            background=Palette.ACCENT_DARK,
            foreground="white",
            relief="flat",
            font=("Segoe UI", 10, "bold"),
        )
        style.map(
            "AnalogSan.Treeview",
            background=[("selected", Palette.SEL)],
            foreground=[("selected", "white")],
        )
        style.configure("Vertical.TScrollbar",   background=Palette.BG_CARD,
                         troughcolor=Palette.BG_CANVAS, borderwidth=0)
        style.configure("Horizontal.TScrollbar", background=Palette.BG_CARD,
                         troughcolor=Palette.BG_CANVAS, borderwidth=0)

    def _canvas_hint(self):
        """Exibe mensagem inicial no canvas de visualização."""
        self._mol_canvas.delete("all")
        self._mol_canvas.create_text(
            158, 105,
            text="Cole um SMILES acima\npara visualizar a molécula",
            fill=Palette.TEXT_DIM, font=("Segoe UI", 11),
            justify="center",
        )

    # ---------------------------------------------------------------------
    #  CALLBACKS - EVENTOS DE UI
    # ---------------------------------------------------------------------

    def _toggle_theme(self):
        if ctk.get_appearance_mode() == "Dark":
            ctk.set_appearance_mode("light")
            self._theme_btn.configure(text="Modo Escuro")
        else:
            ctk.set_appearance_mode("dark")
            self._theme_btn.configure(text="Modo Claro")

    def _load_example(self, smiles: str):
        self._smiles_entry.delete(0, "end")
        self._smiles_entry.insert(0, smiles)
        self._render_molecule_canvas(smiles, self._mol_canvas, size=(316, 210))

    def _on_smiles_change(self, event=None):
        smiles = self._smiles_entry.get().strip()
        if smiles:
            self._render_molecule_canvas(smiles, self._mol_canvas, size=(316, 210))
        else:
            self._canvas_hint()

    def _on_qty_slider(self, value):
        qty = int(float(value))
        self._qty_val_lbl.configure(text=str(qty))
        self._qty_entry.delete(0, "end")
        self._qty_entry.insert(0, str(qty))

    def _on_qty_entry_enter(self, event=None):
        try:
            val = int(self._qty_entry.get().strip())
            val = max(1, min(1000, val))
        except ValueError:
            val = 100
        self._qty_slider.set(val)
        self._qty_val_lbl.configure(text=str(val))
        self._qty_entry.delete(0, "end")
        self._qty_entry.insert(0, str(val))

    def _on_row_select(self, event=None):
        """Exibe detalhes e miniatura da molécula selecionada na tabela."""
        sel = self._table.selection()
        if not sel:
            return
        item = self._table.item(sel[0])
        vals = item["values"]
        if not vals:
            return

        analog_id  = vals[0]
        smiles_trunc = vals[1]

        # Recupera o análogo completo pelo ID
        analog = next((a for a in self._analogs if a["ID_Molecula"] == analog_id), None)
        if not analog:
            return

        alerts_text = str(analog.get("AlertasQuimicos", "n/d"))
        if len(alerts_text) > 92:
            alerts_text = alerts_text[:89] + "..."
        evidence_type = {
            "reaction_class_hypothesis": "classe de reação hipotética",
            "graph_edit": "edição de grafo; rota não codificada",
            "scaffold_hop": "troca de scaffold; rota não codificada",
        }.get(analog.get("EvidenciaTipoRegistro"), "n/d")

        objective_detail = (
            f"Frente Pareto: {analog.get('FrentePareto')}\n"
            f"Objetivos    : {analog.get('Objetivos')}\n"
            if "FrentePareto" in analog else ""
        )
        detail = (
            f"ID         : {analog['ID_Molecula']}\n"
            f"Modificação: {analog['Modificacao']}\n"
            f"Prior. regra: R{analog.get('PrioridadeRegra', '')} "
            f"(nível manual: {analog.get('NivelRegraManual', 'n/d')})\n"
            f"Classe (hip.): {analog.get('EvidenciaClasseReacao', 'n/d')}\n"
            f"Tipo de reg.: {evidence_type}\n"
            f"Evidência   : {'não catalogada' if analog.get('EvidenciaStatus') == 'not_catalogued' else analog.get('EvidenciaStatus', 'n/d')} "
            f"({analog.get('EvidenciaFonteCount', 0)} fontes; rota verificada: {'sim' if analog.get('EvidenciaRotaVerificada') else 'não'})\n"
            f"Nota        : {analog.get('EvidenciaObservacao', 'n/d')}\n"
            f"{objective_detail}"
            f"Origem     : {analog.get('EtapaGeracao', 'n/d')}\n"
            f"Score heur.: {analog.get('ScoreGlobal', 'n/d')} "
            f"({analog.get('CategoriaScore', 'n/d')})\n"
            f"Tanimoto   : {analog.get('SimilaridadeTanimoto', 'n/d')}\n"
            f"SA Score   : {analog.get('SAScore', 'n/d')}\n"
            f"Alertas    : {analog.get('QtdAlertas', 'n/d')} "
            f"({alerts_text})\n"
            f"MW         : {analog['MW']} Da\n"
            f"LogP       : {analog['LogP']}\n"
            f"HBA / HBD  : {analog['HBA']} / {analog['HBD']}\n"
            f"TPSA       : {analog['TPSA']} A^2\n"
            f"RotBonds   : {analog['RotBonds']}\n"
            f"HAC        : {analog['HAC']}\n"
            f"\nSMILES:\n{analog['SMILES']}"
        )
        self._detail_lbl.configure(text=detail, text_color=Palette.TEXT)
        self._render_molecule_canvas(
            analog["SMILES"], self._detail_canvas, size=(220, 130)
        )

    # ---------------------------------------------------------------------
    #  CALLBACKS - AÇÕES PRINCIPAIS
    # ---------------------------------------------------------------------

    def _selected_routes(self) -> List[str]:
        routes = []
        if self._classic_route_var.get():
            routes.append(ce.ROUTE_CLASSIC)
        if self._non_classic_route_var.get():
            routes.append(ce.ROUTE_NON_CLASSIC)
        return routes

    def _selected_objectives(self) -> List[Dict[str, str]]:
        direction_map = {"Aumentar": "increase", "Reduzir": "decrease", "Manter": "maintain"}
        return [
            {"property": prop.get(), "direction": direction_map[direction.get()]}
            for prop, direction in self._objective_vars if prop.get() != "Nenhum"
        ]

    def _on_generate(self):
        smiles = self._smiles_entry.get().strip()
        if not smiles:
            messagebox.showwarning("Atenção", "Cole um código SMILES válido no campo de entrada.")
            return
        if not ce.validate_smiles(smiles):
            messagebox.showerror("SMILES Inválido",
                                  f"Não foi possível interpretar o SMILES:\n\n{smiles}\n\n"
                                  "Verifique a sintaxe e tente novamente.")
            return

        self._on_qty_entry_enter()
        try:
            qty = int(float(self._qty_slider.get()))
        except Exception:
            qty = 100

        selected_routes = self._selected_routes()
        if not selected_routes:
            messagebox.showwarning(
                "Atenção",
                "Selecione pelo menos uma rota de modificação molecular.",
            )
            return

        objective_profile = self._selected_objectives()
        selected_properties = [item["property"] for item in objective_profile]
        if len(selected_properties) != len(set(selected_properties)):
            messagebox.showwarning("Objetivo repetido", "Escolha propriedades diferentes nos objetivos.")
            return

        try:
            max_sa_score = float(self._max_sa_var.get().strip())
            if not 1.0 <= max_sa_score <= 10.0:
                raise ValueError
        except ValueError:
            messagebox.showwarning("Parâmetro inválido", "O SA Score máximo deve ser um número entre 1 e 10.")
            return

        # Reseta estado
        self._clear_table()
        self._clear_log()
        self._analogs = []
        self._stop_requested = False
        self._running = True

        self._gen_btn.configure(state="disabled")
        self._stop_btn.configure(state="normal")
        self._export_btn.configure(state="disabled")
        self._report_btn.configure(state="disabled")
        self._progress_bar.set(0)
        self._progress_lbl.configure(
            text=f"Iniciando geração de {qty} análogos...",
            text_color=Palette.TEXT_DIM,
        )
        self._stats_lbl.configure(text="")
        self._last_base_smiles = smiles
        self._last_generation_params = {
            "target_count": qty,
            "use_veber": self._veber_var.get(),
            "selected_routes": selected_routes,
            "max_sa_score": max_sa_score,
            "allow_multistep": self._multistep_var.get(),
            "seed": 0xA11A,
            "objective_profile": objective_profile,
        }

        self._log(f"--- Sessão iniciada - {qty} análogos solicitados ---")

        # Inicia thread de geração (daemon: encerra com a janela)
        t = threading.Thread(
            target=self._generation_worker,
            args=(smiles, qty, self._veber_var.get(), selected_routes, max_sa_score, self._multistep_var.get()),
            daemon=True,
        )
        t.start()

    def _on_stop(self):
        self._stop_requested = True
        self._log("Sinal de parada enviado - aguardando encerramento...")
        self._stop_btn.configure(state="disabled")

    def _on_export(self):
        if not self._analogs:
            messagebox.showwarning("Atenção", "Nenhum análogo disponível para exportação.")
            return

        filepath = filedialog.asksaveasfilename(
            title="Salvar Análogos",
            defaultextension=".csv",
            filetypes=[
                ("CSV - Docking / ADMET (recomendado)", "*.csv"),
                ("SMILES - Virtual Screening",          "*.smi"),
                ("MOL2 - Coordenadas 3D / docking",     "*.mol2"),
                ("Todos os arquivos",                    "*.*"),
            ],
            initialfile="analog_san_analogs",
        )
        if not filepath:
            return

        try:
            if filepath.lower().endswith(".smi"):
                ce.export_to_smi(self._analogs, filepath)
                fmt = "SMILES (.smi)"
            elif filepath.lower().endswith(".mol2"):
                ce.export_to_mol2(self._analogs, filepath)
                fmt = "MOL2 (.mol2)"
            else:
                ce.export_to_csv(self._analogs, filepath)
                fmt = "CSV"

            self._log(f"Exportado ({fmt}): {filepath}")
            messagebox.showinfo(
                "Exportação Concluída",
                f"{len(self._analogs)} análogos exportados com sucesso!\n\n"
                f"Arquivo: {filepath}\n\n"
                "Dica: arquivos MOL2 já incluem coordenadas 3D para preparo\n"
                "posterior em fluxos de docking molecular.",
            )
        except Exception as e:
            messagebox.showerror("Erro na Exportação", str(e))

    def _on_report(self):
        if not self._analogs:
            messagebox.showwarning("Atenção", "Nenhum análogo disponível para relatório.")
            return

        output_dir = filedialog.askdirectory(
            title="Selecionar pasta para relatório científico",
        )
        if not output_dir:
            return

        try:
            manifest = ce.export_validation_report(
                analogs=self._analogs,
                base_smiles=self._last_base_smiles or self._smiles_entry.get().strip(),
                output_dir=output_dir,
                params=self._last_generation_params,
            )
            self._log(f"Relatório científico gerado: {manifest['output_dir']}")
            self._log(f"Hash dos resultados: {manifest['result_sha256']}")
            messagebox.showinfo(
                "Relatório Gerado",
                "Relatório científico, CSV validado, metadados e figuras PNG/SVG "
                "foram gerados com sucesso.\n\n"
                f"Pasta: {manifest['output_dir']}\n"
                f"Melhor candidato: {manifest['top_id']} "
                f"(Score {manifest['top_score']})",
            )
        except Exception as e:
            messagebox.showerror("Erro no Relatório", str(e))

    # ---------------------------------------------------------------------
    #  WORKER THREAD E QUEUE POLLING
    # ---------------------------------------------------------------------

    def _generation_worker(
        self,
        smiles: str,
        qty: int,
        use_veber: bool,
        selected_routes: List[str],
        max_sa_score: float,
        allow_multistep: bool,
    ):
        """
        Executa a geração em thread separada.
        Comunica-se com a GUI exclusivamente via self._gui_queue.
        """
        try:
            analogs = ce.generate_analogs(
                base_smiles=smiles,
                target_count=qty,
                use_veber=use_veber,
                selected_routes=selected_routes,
                max_sa_score=max_sa_score,
                allow_multistep=allow_multistep,
                progress_callback=lambda cur, tot:
                    self._gui_queue.put(("progress", (cur, tot))),
                log_callback=lambda msg:
                    self._gui_queue.put(("log", msg)),
                stop_flag=lambda: self._stop_requested,
                random_seed=self._last_generation_params.get("seed", 0xA11A),
                objective_profile=self._last_generation_params.get("objective_profile"),
            )
            self._gui_queue.put(("done", analogs))
        except Exception as exc:
            self._gui_queue.put(("error", str(exc)))

    def _poll_queue(self):
        """
        Polling da queue a cada 80 ms - única forma thread-safe de atualizar
        widgets tkinter a partir de threads externas.
        """
        try:
            while True:
                msg_type, payload = self._gui_queue.get_nowait()
                if msg_type == "progress":
                    self._handle_progress(payload)
                elif msg_type == "log":
                    self._log(payload)
                elif msg_type == "done":
                    self._handle_done(payload)
                elif msg_type == "error":
                    self._handle_error(payload)
        except queue.Empty:
            pass
        finally:
            self.after(80, self._poll_queue)

    def _handle_progress(self, payload):
        cur, tot = payload
        pct = cur / tot if tot > 0 else 0
        self._progress_bar.set(pct)
        self._progress_lbl.configure(
            text=f"{cur}/{tot} aprovados  ({pct*100:.0f}%)",
            text_color=Palette.TEXT_GREEN if pct > 0 else Palette.TEXT_DIM,
        )

    def _handle_done(self, analogs: List[Dict]):
        self._analogs = analogs
        self._running = False
        self._populate_table(analogs)
        requested = int(self._last_generation_params.get("target_count", len(analogs)))
        self._progress_bar.set(min(1.0, len(analogs) / requested) if requested else 1.0)
        self._progress_lbl.configure(
            text=(f"{len(analogs)}/{requested} candidatos aprovados" if len(analogs) < requested
                  else f"{len(analogs)} candidatos aprovados"),
            text_color=Palette.TEXT_GREEN if len(analogs) >= requested else Palette.TEXT_YELLOW,
        )
        self._stats_lbl.configure(
            text=f"Total: {len(analogs)} mol  |  "
                 f"Score médio: {sum(a.get('ScoreGlobal', 0) for a in analogs)/len(analogs):.1f}  |  "
                 f"LogP médio: {sum(a['LogP'] for a in analogs)/len(analogs):.2f}  |  "
                 f"MW médio: {sum(a['MW'] for a in analogs)/len(analogs):.1f} Da"
            if analogs else ""
        )
        self._gen_btn.configure(state="normal")
        self._stop_btn.configure(state="disabled")
        if analogs:
            self._export_btn.configure(state="normal")
            self._report_btn.configure(state="normal")
        messagebox.showinfo(
            "Geração Concluída",
            f"{len(analogs)} candidatos passaram pelos filtros configurados (solicitados: {requested}).\n"
            f"A viabilidade de rota sintética não foi avaliada pelo programa.\n\n"
            f"Use 'Exportar' para salvar CSV/SMI/MOL2 ou 'Relatório' para gerar figuras científicas.",
        )

    def _handle_error(self, msg: str):
        self._running = False
        self._gen_btn.configure(state="normal")
        self._stop_btn.configure(state="disabled")
        self._report_btn.configure(state="disabled")
        self._log(f"ERRO: {msg}")
        messagebox.showerror("Erro na Geração", f"Ocorreu um erro durante a geração:\n\n{msg}")

    # ---------------------------------------------------------------------
    #  MANIPULAÇÃO DA TABELA
    # ---------------------------------------------------------------------

    def _populate_table(self, analogs: List[Dict]):
        """Preenche a Treeview com todos os análogos aprovados."""
        self._clear_table()
        for i, a in enumerate(analogs):
            tag = "even" if i % 2 == 0 else "odd"
            smiles_display = a["SMILES"] if len(a["SMILES"]) <= 55 else a["SMILES"][:52] + "..."
            self._table.insert("", "end", iid=a["ID_Molecula"], values=(
                a["ID_Molecula"],
                smiles_display,
                a["Modificacao"],
                f"R{a.get('PrioridadeRegra', '')} {a.get('NivelRegraManual', '')}",
                a.get("FrentePareto", ""),
                a.get("Objetivos", ""),
                f"{a.get('ScoreGlobal', 0):.1f}",
                f"{a.get('SimilaridadeTanimoto', 0):.3f}",
                f"{a.get('SAScore', 0):.2f}",
                a.get("QtdAlertas", 0),
                f"{a['MW']:.1f}",
                f"{a['LogP']:.2f}",
                a["HBA"],
                a["HBD"],
                f"{a['TPSA']:.1f}",
                a["RotBonds"],
            ), tags=(tag,))

    def _clear_table(self):
        for item in self._table.get_children():
            self._table.delete(item)

    def _on_table_click(self, event):
        """Ordena ao clicar no cabeçalho sem usar a opção command do ttk."""
        if self._table.identify_region(event.x, event.y) != "heading":
            return None

        column_id = self._table.identify_column(event.x)
        if not column_id.startswith("#"):
            return None

        try:
            column_index = int(column_id[1:]) - 1
            column = self._table["columns"][column_index]
        except (ValueError, IndexError):
            return None

        if column in self._sortable_columns:
            self._sort_column(column)
            return "break"

        return None

    def _sort_column(self, col: str):
        """Ordena a tabela pela coluna clicada (alterna asc <-> desc)."""
        rows = [(self._table.set(k, col), k) for k in self._table.get_children("")]
        reverse = self._sort_reverse.get(col, False)
        try:
            rows.sort(key=lambda t: float(t[0].replace("...", "")), reverse=reverse)
        except ValueError:
            rows.sort(key=lambda t: t[0], reverse=reverse)
        for idx, (_, k) in enumerate(rows):
            self._table.move(k, "", idx)
            self._table.item(k, tags=("even" if idx % 2 == 0 else "odd",))
        self._sort_reverse[col] = not reverse

    # ---------------------------------------------------------------------
    #  VISUALIZAÇÃO DE MOLÉCULAS
    # ---------------------------------------------------------------------

    def _render_molecule_canvas(
        self,
        smiles: str,
        canvas: tk.Canvas,
        size: tuple = (316, 210),
    ):
        """
        Renderiza um SMILES no canvas tk especificado.
        Usa chem_engine.mol_to_pil_image() para geração da imagem.
        """
        img = ce.mol_to_pil_image(smiles, size=size)
        if img:
            tk_img = ImageTk.PhotoImage(img)
            canvas.delete("all")
            canvas.create_image(size[0] // 2, size[1] // 2, image=tk_img)
            # Armazena referência para evitar coleta pelo GC do Python
            if canvas is self._mol_canvas:
                self._canvas_img_ref = tk_img
            else:
                self._detail_img_ref = tk_img
        else:
            canvas.delete("all")
            canvas.create_text(
                size[0] // 2, size[1] // 2,
                text="SMILES inválido",
                fill=Palette.ACCENT_RED, font=("Segoe UI", 11),
            )

    # ---------------------------------------------------------------------
    #  LOG
    # ---------------------------------------------------------------------

    def _log(self, message: str):
        self._log_box.configure(state="normal")
        self._log_box.insert("end", message + "\n")
        self._log_box.see("end")
        self._log_box.configure(state="disabled")

    def _clear_log(self):
        self._log_box.configure(state="normal")
        self._log_box.delete("1.0", "end")
        self._log_box.configure(state="disabled")
