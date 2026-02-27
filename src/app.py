import streamlit as st
import tempfile
import os
import io
import pandas as pd
import plotly.graph_objects as go

# Importiamo le funzioni dal file sorgente
from src_armSol import carica_input_excel, crea_armatura_definitiva, crea_geometria_impalcato

# ==========================================
# 1. MOTORE GRAFICO PLOTLY (PONTE + TABELLA INTEGRATA)
# ==========================================

def crea_grafico_ponte_integrato(conci_estratti, campate_mm):
    """Genera il disegno del ponte con la tabella dei dati perfettamente allineata sotto di esso."""
    fig = go.Figure()

    lung_assi = sum(campate_mm)
    max_h = max([c['h_start'] for c in conci_estratti] + [c['h_end'] for c in conci_estratti])

    # ---------------------------------------------------------
    # PARTE 1: DISEGNO DEL PONTE
    # ---------------------------------------------------------
    for c in conci_estratti:
        x_coords = [c['x1'], c['x2'], c['x2'], c['x1'], c['x1']]
        y_coords = [0, 0, -c['h_end'], -c['h_start'], 0]
        
        hover_testo = (
            f"<b>{c['nome']}</b><br>Lunghezza: {c['L_mm']/1000} m<br>"
            f"H Iniziale: {c['h_start']} mm | H Finale: {c['h_end']} mm<br>"
            f"Armatura Sup: Φ{c['phi_sup']}/{c['passo_sup_cm']} cm"
        )
        if c['phi_inf']: hover_testo += f"<br>Armatura Inf: Φ{c['phi_inf']}/{c['passo_inf_cm']} cm"

        fig.add_trace(go.Scatter(
            x=x_coords, y=y_coords, fill='toself', fillcolor='rgba(169, 169, 169, 0.4)', 
            line=dict(color='black', width=1.5), name=c['nome'], hoverinfo='text', text=hover_testo, showlegend=False
        ))

        fig.add_trace(go.Scatter(
            x=[(c['x1'] + c['x2']) / 2], y=[(-c['h_start'] - c['h_end']) / 4], 
            mode='text', text=[c['nome']], textfont=dict(color='black', size=11), hoverinfo='skip', showlegend=False
        ))

    supporti_x = [0]
    x_curr = 0
    for L in campate_mm: x_curr += L; supporti_x.append(x_curr)
    
    for idx, sx in enumerate(supporti_x):
        nome_appoggio = "SPA" if idx == 0 else "SPB" if idx == len(supporti_x) - 1 else f"P{idx}"
        fig.add_trace(go.Scatter(x=[sx, sx], y=[500, -max_h - 2000], mode='lines', line=dict(color='red', width=2, dash='dashdot'), hoverinfo='skip', showlegend=False))
        fig.add_trace(go.Scatter(x=[sx], y=[1500], mode='text', text=[f"<b>{nome_appoggio}</b>"], textfont=dict(color='red', size=13), hoverinfo='skip', showlegend=False))

    fig.add_trace(go.Scatter(x=[0, lung_assi], y=[0, 0], mode='lines', line=dict(color='blue', width=1.5, dash='dash'), hoverinfo='skip', showlegend=False))


    # ---------------------------------------------------------
    # PARTE 2: TABELLA INTEGRATA (Stile CAD senza prima colonna)
    # ---------------------------------------------------------
    y_tab_start = -max_h - 4000  
    h_row = 1500  
    w_header = 25000  # Spazio riservato a sinistra per i titoli fluttuanti
    
    righe = [
        ("Concio Geometrico", "nome"),
        ("Tipo concio strutturale", "strutturale"),
        ("Lunghezza concio (mm)", "L_mm"),
        ("Pioli Nelson", "pioli"),
        ("Piattabanda sup.", "psup"),
        ("Saldatura anima-pb.sup.", "sald_psup"),
        ("H. iniziale (mm)", "h_start"),
        ("H. finale (mm)", "h_end"),
        ("Spessore anima", "sa"),
        ("Saldatura anima-pb.inf.", "sald_pinf"),
        ("Piattabanda inf.", "pinf"),
        ("Armatura Sup. Soletta", "arm_sup"),
        ("Armatura Inf. Soletta", "arm_inf")
    ]
    
    y_tab_end = y_tab_start - len(righe) * h_row

    # --- Disegno Griglia Tabella (Linee) ---
    shapes = []
    
    # Tutte le linee Orizzontali partono ora da X=0 (non più da -w_header)
    for i in range(len(righe) + 1):
        y_line = y_tab_start - i * h_row
        width = 2 if i == 0 or i == len(righe) else 0.5
        shapes.append(dict(type="line", x0=0, y0=y_line, x1=lung_assi, y1=y_line, line=dict(color="black", width=width)))
        
    # Linee Verticali (ora solo quella a X=0 e la chiusura finale)
    shapes.append(dict(type="line", x0=0, y0=y_tab_start, x1=0, y1=y_tab_end, line=dict(color="black", width=2)))
    shapes.append(dict(type="line", x0=lung_assi, y0=y_tab_start, x1=lung_assi, y1=y_tab_end, line=dict(color="black", width=2)))
    
    # Linee Verticali interne (Allineate ai Conci)
    for c in conci_estratti:
        if c['x2'] < lung_assi:
            shapes.append(dict(type="line", x0=c['x2'], y0=y_tab_start, x1=c['x2'], y1=y_tab_end, line=dict(color="gray", width=0.5, dash='dot')))

    fig.update_layout(shapes=shapes)

    # --- Compilazione Testi Tabella ---
    titoli_x, titoli_y, titoli_text = [], [], []
    celle_x, celle_y, celle_text = [], [], []

    for i, r in enumerate(righe):
        y_text = y_tab_start - i * h_row - (h_row / 2)
        
        # Testo Colonna Titoli: si posizionano poco prima di X=0 e si estendono verso sinistra ("middle left" significa punto a destra, testo a sinistra)
        titoli_x.append(-1000) 
        titoli_y.append(y_text)
        titoli_text.append(f"<b>{r[0]}</b>")
        
        for c in conci_estratti:
            x_mid = (c['x1'] + c['x2']) / 2
            
            if r[1] == 'arm_sup': val = f"Φ{c.get('phi_sup', '-')}/{c.get('passo_sup_cm', '-')} cm"
            elif r[1] == 'arm_inf': val = f"Φ{c.get('phi_inf', '-')}/{c.get('passo_inf_cm', '-')} cm" if c.get('phi_inf') else "-"
            else: 
                v = str(c.get(r[1], '-')).replace('.0', '')
                val = v if v not in ['nan', 'None', ''] else "-"
                
            celle_x.append(x_mid)
            celle_y.append(y_text)
            celle_text.append(val)

    # Inserimento Titoli (Fluttuanti)
    fig.add_trace(go.Scatter(
        x=titoli_x, y=titoli_y, mode='text', text=titoli_text,
        textposition="middle left", textfont=dict(color='black', size=11), hoverinfo='skip', showlegend=False
    ))
    
    # Inserimento Celle
    fig.add_trace(go.Scatter(
        x=celle_x, y=celle_y, mode='text', text=celle_text,
        textposition="middle center", textfont=dict(color='black', size=11), hoverinfo='skip', showlegend=False
    ))

    # Punto invisibile per costringere Plotly a mostrare i testi fluttuanti a sinistra
    fig.add_trace(go.Scatter(x=[-w_header], y=[y_tab_end], mode='markers', marker=dict(color='rgba(0,0,0,0)'), hoverinfo='skip', showlegend=False))

    fig.update_layout(
        title="Vista Impalcato e Tabella Dati Sincronizzata (Usa la rotellina per lo Zoom)",
        plot_bgcolor='white', margin=dict(l=10, r=10, t=40, b=10),
        xaxis=dict(showgrid=False, zeroline=False, visible=False),
        yaxis=dict(showgrid=False, zeroline=False, visible=False),
        height=750, 
        hovermode="closest"
    )
    
    return fig


# ==========================================
# 2. CONFIGURAZIONE PAGINA STREAMLIT
# ==========================================
st.set_page_config(page_title="Ponti - Generatore DXF", page_icon="🌉", layout="wide")

st.title("🌉 Generatore Elaborati Impalcato e Armature")
st.markdown("Carica il file Excel dei conci. L'anteprima CAD apparirà **istantaneamente**. Clicca su *Genera* per esportare.")

for key in ['geometria_bytes', 'armature_bytes', 'tabella_bytes']:
    if key not in st.session_state: st.session_state[key] = None

# --- UPLOAD DEL FILE ---
uploaded_file = st.file_uploader("Trascina qui il file Excel (.xlsx)", type=['xlsx'])

if uploaded_file is not None:
    with tempfile.TemporaryDirectory() as tmpdir:
        input_path = os.path.join(tmpdir, "input.xlsx")
        with open(input_path, "wb") as f: f.write(uploaded_file.getbuffer())

        # ===============================================
        # 3. LETTURA ED ANTEPRIMA IMMEDIATA
        # ===============================================
        conci_estratti, _, _, campate_lette = carica_input_excel(input_path)

        if not campate_lette:
            st.error("❌ ERRORE: Assicurati che esista il foglio 'Campate' nel tuo file Excel.")
        else:
            st.success("✅ File letto con successo! Visualizza l'anteprima in basso.")
            
            # Mostra il super-grafico integrato
            st.plotly_chart(crea_grafico_ponte_integrato(conci_estratti, campate_lette), use_container_width=True)

            # ===============================================
            # 4. BOTTONE DI ESPORTAZIONE
            # ===============================================
            st.markdown("---")
            st.subheader("💾 Esportazione File CAD ed Excel")
            
            if st.button("🚀 Genera DXF e Distinte", type="primary"):
                with st.spinner("Calcolo vettoriale e Nesting in corso..."):
                    out_armature = os.path.join(tmpdir, "ponte_armature.dxf")
                    out_geometria = os.path.join(tmpdir, "ponte_geometria.dxf")
                    out_tabella = os.path.join(tmpdir, "distinta_e_ordini.xlsx")

                    try:
                        crea_armatura_definitiva(out_armature, input_path, out_tabella, campate_lette)
                        crea_geometria_impalcato(out_geometria, input_path, campate_lette)
                        
                        with open(out_geometria, "rb") as f: st.session_state.geometria_bytes = f.read()
                        with open(out_armature, "rb") as f: st.session_state.armature_bytes = f.read()
                        with open(out_tabella, "rb") as f: st.session_state.tabella_bytes = f.read()
                    except Exception as e:
                        st.error(f"❌ Errore durante la generazione: {e}")

            # ===============================================
            # 5. AREA DOWNLOAD
            # ===============================================
            if st.session_state.geometria_bytes:
                st.success("✅ Generazione completata! Scarica i file dai pulsanti sottostanti.")
                
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.download_button("📐 Scarica DXF Geometria", data=st.session_state.geometria_bytes, file_name="ponte_geometria.dxf", mime="application/dxf", use_container_width=True)
                with col2:
                    st.download_button("🏗️ Scarica DXF Armature", data=st.session_state.armature_bytes, file_name="ponte_armature.dxf", mime="application/dxf", use_container_width=True)
                with col3:
                    st.download_button("📊 Scarica Excel Distinta", data=st.session_state.tabella_bytes, file_name="distinta_e_ordini.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", use_container_width=True)