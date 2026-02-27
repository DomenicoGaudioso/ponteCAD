import streamlit as st
import tempfile
import os

# Importiamo le funzioni dal file sorgente
from src_armSol import carica_input_excel, crea_armatura_definitiva, crea_geometria_impalcato

# --- CONFIGURAZIONE PAGINA ---
st.set_page_config(
    page_title="Ponti - DXF",
    page_icon="🌉",
    layout="centered"
)

st.title("🌉 Prospetto Impalcato e Armature")
st.markdown("""
Questo strumento genera automaticamente **Tavole Costruttive DXF** e le **Distinte Ferri in Excel** partendo da un foglio di calcolo. 
Assicurati che il tuo file Excel contenga i fogli **'Conci'** e **'Campate'**.
""")

# Inizializziamo lo stato della sessione per conservare i file in memoria
if 'geometria_bytes' not in st.session_state:
    st.session_state.geometria_bytes = None
if 'armature_bytes' not in st.session_state:
    st.session_state.armature_bytes = None
if 'tabella_bytes' not in st.session_state:
    st.session_state.tabella_bytes = None
if 'generazione_completata' not in st.session_state:
    st.session_state.generazione_completata = False

# --- UPLOAD DEL FILE ---
uploaded_file = st.file_uploader("Trascina qui il tuo file Excel di Input (.xlsx)", type=['xlsx'])

if uploaded_file is not None:
    st.success("File caricato correttamente!")
    
    # Bottone di avvio
    if st.button("🚀 Genera Elaborati", type="primary"):
        with st.spinner("Elaborazione in corso, attendere..."):
            
            # Creiamo una cartella temporanea per manipolare i file
            with tempfile.TemporaryDirectory() as tmpdir:
                
                input_path = os.path.join(tmpdir, "input.xlsx")
                with open(input_path, "wb") as f:
                    f.write(uploaded_file.getbuffer())

                # 1. Lettura propedeutica per controllo validità file
                _, _, _, campate_lette = carica_input_excel(input_path)

                if not campate_lette:
                    st.error("❌ ERRORE: Non è stato possibile leggere le campate. Assicurati che esista il foglio 'Campate'.")
                    st.session_state.generazione_completata = False
                else:
                    out_armature = os.path.join(tmpdir, "ponte_armature.dxf")
                    out_geometria = os.path.join(tmpdir, "ponte_geometria.dxf")
                    out_tabella = os.path.join(tmpdir, "distinta_e_ordini.xlsx")

                    # 2. ESECUZIONE ALGORITMI
                    try:
                        crea_armatura_definitiva(out_armature, input_path, out_tabella, campate_lette)
                        crea_geometria_impalcato(out_geometria, input_path, campate_lette)
                        
                        # 3. SALVATAGGIO DEI FILE NELLA MEMORIA DI STREAMLIT
                        with open(out_geometria, "rb") as f:
                            st.session_state.geometria_bytes = f.read()
                        with open(out_armature, "rb") as f:
                            st.session_state.armature_bytes = f.read()
                        with open(out_tabella, "rb") as f:
                            st.session_state.tabella_bytes = f.read()
                            
                        # Diciamo al programma che abbiamo finito con successo
                        st.session_state.generazione_completata = True
                            
                    except Exception as e:
                        st.error(f"❌ Si è verificato un errore durante la generazione: {e}")
                        st.session_state.generazione_completata = False

    # --- BOX DI DOWNLOAD (Slegato dal pulsante "Genera") ---
    # Mostra i pulsanti solo se i file sono presenti in memoria
    if st.session_state.generazione_completata:
        st.success("✅ Generazione completata con successo! I file sono pronti per il download.")
        st.markdown("### 📥 Download File Generati")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.download_button(
                label="📐 DXF Geometria",
                data=st.session_state.geometria_bytes,
                file_name="ponte_geometria.dxf",
                mime="application/dxf"
            )
            
        with col2:
            st.download_button(
                label="🏗️ DXF Armature",
                data=st.session_state.armature_bytes,
                file_name="ponte_armature.dxf",
                mime="application/dxf"
            )
            
        with col3:
            st.download_button(
                label="📊 Excel Distinta",
                data=st.session_state.tabella_bytes,
                file_name="distinta_e_ordini.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )