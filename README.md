# ponteCAD
# 🌉 Generatore Elaborati e Armature per Ponti (Bridge CAD Generator)

Un'applicazione web sviluppata in **Python** e **Streamlit** per l'automazione della progettazione strutturale di ponti a sezione mista o in CAP. 
Il software legge un foglio di calcolo Excel contenente la parametrizzazione geometrica dell'impalcato e genera istantaneamente **tavole costruttive in formato DXF** e **computi metrici avanzati in Excel**.

---

## ✨ Funzionalità Principali

* **📐 DXF Geometria:** Generazione del tracciato longitudinale del ponte, con quotatura automatica progressiva (conci, campate, sviluppo totale). Disegna in automatico il variare delle altezze dell'impalcato e impagina una tabella riassuntiva dinamica sotto l'asse.
* **🏗️ DXF Armature (Soletta):** Disegno automatico delle barre di armatura longitudinali superiori e inferiori. Calcola in automatico:
  * Sovrapposizioni di normativa (min. 50 diametri).
  * Sfalsamento dei giunti (Staggering) tra l'armatura superiore e inferiore.
  * Fasi di getto della soletta con simbologia CAD.
* **📊 Distinta Ferri e Computo:** Esportazione in Excel di tutti i ferri calcolati con indicazione di posizione, diametro, lunghezza, numero di pezzi e peso complessivo.
* **✂️ Ottimizzazione Sfridi (Nesting FFD):** Algoritmo integrato (First-Fit Decreasing) che calcola esattamente quante barre commerciali da 12 metri ordinare per ridurre al minimo gli sfridi di cantiere.
* **🌐 Interfaccia Web:** Dashboard intuitiva basata su Streamlit per caricare i file e scaricare i risultati in un clic, rendendo il tool fruibile anche a chi non sa programmare.

---

## 📂 Struttura del Progetto

```text
├── app.py                # Frontend: Interfaccia utente web creata con Streamlit
├── src_armSol.py         # Backend: Motore di calcolo e generazione file DXF/Excel (ezdxf, pandas)
├── requirements.txt      # Dipendenze e librerie del progetto
└── input_conci.xlsx      # Esempio di file Excel di input (Template)
