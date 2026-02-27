import ezdxf
from ezdxf.enums import TextEntityAlignment
import pandas as pd
import math

# ==========================================
# 1. FUNZIONI DI LETTURA E CALCOLO
# ==========================================

def carica_input_excel(filepath):
    # ---------------------------------------------------------
    # 1. LETTURA FOGLIO CONCI (Cerca foglio 'Conci' o prende il 1°)
    # ---------------------------------------------------------
    try: df_conci = pd.read_excel(filepath, sheet_name='Conci')
    except ValueError: df_conci = pd.read_excel(filepath, sheet_name=0)
    
    conci = []
    x_corr = 0.0
    disegna_inf = False
    
    def safe_str(row_data, col_name):
        if col_name in df_conci.columns and pd.notna(row_data[col_name]):
            val = str(row_data[col_name]).strip()
            return val if val != "" else "-"
        return "-"

    for _, row in df_conci.iterrows():
        L_mm = float(row['Lunghezza']) * 1000.0 if 'Lunghezza' in df_conci.columns else 12000.0
        
        try: phi_sup = int(float(row['Φsup'])) if 'Φsup' in df_conci.columns else 20
        except: phi_sup = 20
        try: passo_sup = int(float(row['passo_sup'])) if 'passo_sup' in df_conci.columns else 200
        except: passo_sup = 200

        phi_inf, passo_inf = None, None
        if 'Φinf' in df_conci.columns:
            try:
                val_phi = str(row['Φinf']).strip().lower()
                if val_phi not in ['none', 'nan', '']: phi_inf = int(float(val_phi))
            except: pass
            
        if 'passo_inf' in df_conci.columns:
            try:
                val_passo = str(row['passo_inf']).strip().lower()
                if val_passo not in ['none', 'nan', '']: passo_inf = int(float(val_passo))
            except: pass

        if phi_inf is not None and passo_inf is not None:
            disegna_inf = True

        try: B_mm = float(row['B']) if 'B' in df_conci.columns else 10000.0
        except: B_mm = 10000.0
        try: h_start = float(row['Hstart']) if 'Hstart' in df_conci.columns else 2000.0 
        except: h_start = 2000.0 
        try: h_end = float(row['Hend']) if 'Hend' in df_conci.columns else 2000.0
        except: h_end = 2000.0

        conci.append({
            'nome': safe_str(row, 'Concio geometrico'),
            'strutturale': safe_str(row, 'Concio Strutturale'),
            'x1': x_corr, 'x2': x_corr + L_mm,
            'L_mm': L_mm,
            'phi_sup': phi_sup, 'passo_sup_cm': int(passo_sup / 10),
            'phi_inf': phi_inf, 'passo_inf_cm': int(passo_inf / 10) if passo_inf else None,
            'B_mm': B_mm, 'h_start': h_start, 'h_end': h_end,
            'psup': safe_str(row, 'psup'), 'sa': safe_str(row, 'sa'), 'pinf': safe_str(row, 'pinf'),
            'pioli': safe_str(row, 'pioli'), 'sald_psup': safe_str(row, 'sald_psup'), 'sald_pinf': safe_str(row, 'sald_pinf')
        })
        x_corr += L_mm
        
    # ---------------------------------------------------------
    # 2. LETTURA FOGLIO CAMPATE (Cerca foglio 'Campate' o prende il 2°)
    # ---------------------------------------------------------
    try:
        df_campate = pd.read_excel(filepath, sheet_name='Campate')
    except ValueError:
        try: df_campate = pd.read_excel(filepath, sheet_name=1)
        except IndexError: 
            print("ERRORE: Non trovo il foglio 'Campate' nell'Excel. Usa un array vuoto.")
            return conci, x_corr, disegna_inf, []

    # Cerca una colonna che assomigli alla lunghezza
    colonna_campate = None
    for col in ['Lunghezza (m)', 'Lunghezza', 'L', 'Campate', 'L(m)']:
        if col in df_campate.columns:
            colonna_campate = col
            break
            
    if colonna_campate:
        campate_mm = [float(val) * 1000.0 for val in df_campate[colonna_campate].dropna()]
    else:
        # Se non trova il nome colonna, prende brutalmente i numeri dalla prima colonna a sinistra
        campate_mm = [float(val) * 1000.0 for val in df_campate.iloc[:, 0].dropna()]

    return conci, x_corr, disegna_inf, campate_mm

def get_prop_da_x(x, conci, strato='sup'):
    """Estrae le proprietà in base alla coordinata. Se l'armatura inf è attiva ma una cella 
    specifica ha 'None', usa un fallback di 20/20 per non far crashare la matematica."""
    def safe_extract(c, s):
        if s == 'sup':
            return {'phi': c['phi_sup'], 'passo_cm': c['passo_sup_cm'], 'B_mm': c['B_mm']}
        else:
            p = c['phi_inf'] if c['phi_inf'] is not None else 20
            pc = c['passo_inf_cm'] if c['passo_inf_cm'] is not None else 20
            return {'phi': p, 'passo_cm': pc, 'B_mm': c['B_mm']}

    for c in conci:
        if c['x1'] - 1e-5 <= x <= c['x2'] + 1e-5:
            return safe_extract(c, strato)
            
    return safe_extract(conci[-1], strato)

def get_sov(phi1, phi2, min_sov_mm, mult):
    return max(min_sov_mm, mult * max(phi1, phi2))

def calcola_posizionamento_barre(campate_mm, conci_estratti, strato='sup', lung_barra_mm=12000, 
                                 min_barra_mm=3000, min_sov_mm=1000, mult_sov=50, shift_sfalsamento_mm=0):
    supporti_x = [0]
    x_corrente = 0
    for L in campate_mm:
        x_corrente += L
        supporti_x.append(x_corrente)
    lung_totale_assi = supporti_x[-1]

    pile_x = supporti_x[1:-1]
    pier_bars = []
    
    for sx in pile_x:
        props = get_prop_da_x(sx, conci_estratti, strato)
        x_start = sx - (lung_barra_mm / 2) + shift_sfalsamento_mm
        x_end = sx + (lung_barra_mm / 2) + shift_sfalsamento_mm
        pier_bars.append({'x1': x_start, 'x2': x_end, 'type': 'pier', 'phi': props['phi']})

    barre_finali = []

    for i in range(len(supporti_x) - 1):
        if i == 0:
            p_bar = pier_bars[0]
            L_sov = get_sov(p_bar['phi'], get_prop_da_x(p_bar['x1'], conci_estratti, strato)['phi'], min_sov_mm, mult_sov)
            X_left, X_right = 0, p_bar['x1'] + L_sov 
            
            curr_r = X_right
            gap_bars = []
            while curr_r - X_left > lung_barra_mm:
                x_start_bar = curr_r - lung_barra_mm
                gap_bars.append({'x1': x_start_bar, 'x2': curr_r})
                phi_c = get_prop_da_x(x_start_bar + lung_barra_mm/2, conci_estratti, strato)['phi']
                phi_n = get_prop_da_x(x_start_bar, conci_estratti, strato)['phi']
                curr_r = x_start_bar + get_sov(phi_c, phi_n, min_sov_mm, mult_sov)

            if curr_r - X_left < min_barra_mm and gap_bars:
                deficit = min_barra_mm - (curr_r - X_left)
                gap_bars[-1]['x1'] += deficit
                curr_r += deficit
            gap_bars.append({'x1': X_left, 'x2': curr_r})

        elif i == len(supporti_x) - 2:
            p_bar = pier_bars[-1]
            L_sov = get_sov(p_bar['phi'], get_prop_da_x(p_bar['x2'], conci_estratti, strato)['phi'], min_sov_mm, mult_sov)
            X_left, X_right = p_bar['x2'] - L_sov, lung_totale_assi
            
            curr_l = X_left
            gap_bars = []
            while X_right - curr_l > lung_barra_mm:
                x_end_bar = curr_l + lung_barra_mm
                gap_bars.append({'x1': curr_l, 'x2': x_end_bar})
                phi_c = get_prop_da_x(curr_l + lung_barra_mm/2, conci_estratti, strato)['phi']
                phi_n = get_prop_da_x(x_end_bar, conci_estratti, strato)['phi']
                curr_l = x_end_bar - get_sov(phi_c, phi_n, min_sov_mm, mult_sov)

            if X_right - curr_l < min_barra_mm and gap_bars:
                deficit = min_barra_mm - (X_right - curr_l)
                gap_bars[-1]['x2'] -= deficit
                curr_l -= deficit
            gap_bars.append({'x1': curr_l, 'x2': X_right})

        else:
            p_left, p_right = pier_bars[i-1], pier_bars[i]
            L_sov_l = get_sov(p_left['phi'], get_prop_da_x(p_left['x2'], conci_estratti, strato)['phi'], min_sov_mm, mult_sov)
            L_sov_r = get_sov(p_right['phi'], get_prop_da_x(p_right['x1'], conci_estratti, strato)['phi'], min_sov_mm, mult_sov)
            
            X_left, X_right = p_left['x2'] - L_sov_l, p_right['x1'] + L_sov_r
            curr_l, curr_r = X_left, X_right
            gap_bars_l, gap_bars_r = [], []

            while curr_r - curr_l > lung_barra_mm:
                x_end_l = curr_l + lung_barra_mm
                gap_bars_l.append({'x1': curr_l, 'x2': x_end_l})
                curr_l = x_end_l - get_sov(get_prop_da_x(curr_l + lung_barra_mm/2, conci_estratti, strato)['phi'], get_prop_da_x(x_end_l, conci_estratti, strato)['phi'], min_sov_mm, mult_sov)
                if curr_r - curr_l <= lung_barra_mm: break

                x_start_r = curr_r - lung_barra_mm
                gap_bars_r.append({'x1': x_start_r, 'x2': curr_r})
                curr_r = x_start_r + get_sov(get_prop_da_x(x_start_r + lung_barra_mm/2, conci_estratti, strato)['phi'], get_prop_da_x(x_start_r, conci_estratti, strato)['phi'], min_sov_mm, mult_sov)

            if curr_r - curr_l < min_barra_mm:
                deficit = min_barra_mm - (curr_r - curr_l)
                if gap_bars_l and gap_bars_r:
                    gap_bars_l[-1]['x2'] -= deficit / 2
                    curr_l -= deficit / 2
                    gap_bars_r[-1]['x1'] += deficit / 2
                    curr_r += deficit / 2
                elif gap_bars_l: gap_bars_l[-1]['x2'] -= deficit; curr_l -= deficit
                elif gap_bars_r: gap_bars_r[-1]['x1'] += deficit; curr_r += deficit
            
            gap_bars = gap_bars_l + [{'x1': curr_l, 'x2': curr_r}] + gap_bars_r

        for gb in gap_bars:
            gb['type'] = 'span'
            props = get_prop_da_x((gb['x1'] + gb['x2'])/2, conci_estratti, strato)
            gb['phi'], gb['passo_cm'], gb['B_mm'] = props['phi'], props['passo_cm'], props['B_mm']
            barre_finali.append(gb)

    for pb in pier_bars:
        props = get_prop_da_x((pb['x1'] + pb['x2'])/2, conci_estratti, strato)
        pb['passo_cm'], pb['B_mm'] = props['passo_cm'], props['B_mm']
        barre_finali.append(pb)

    barre_finali.sort(key=lambda b: b['x1'])
    if barre_finali:
        barre_finali[0]['type'] = 'abutment_left'
        barre_finali[-1]['type'] = 'abutment_right'
        
    return barre_finali

def esegui_nesting_sfridi(dati_tabella, lung_comm_cm=1200):
    riepilogo = []
    diametri_unici = set([r['Diametro (mm)'] for r in dati_tabella if r['Posizione'] != 'TOTALE'])
    
    for phi in sorted(list(diametri_unici)):
        pezzi_phi = [r for r in dati_tabella if r['Diametro (mm)'] == phi and r['Posizione'] != 'TOTALE']
        barre_intere = sum([int(r['N. Barre']) for r in pezzi_phi if r['Lunghezza Barra (cm)'] >= lung_comm_cm])
        
        da_tagliare = []
        for r in pezzi_phi:
            if r['Lunghezza Barra (cm)'] < lung_comm_cm:
                da_tagliare.extend([r['Lunghezza Barra (cm)']] * int(r['N. Barre']))
        
        da_tagliare.sort(reverse=True)
        bins = [] 
        for pezzo in da_tagliare:
            inserito = False
            for i in range(len(bins)):
                if bins[i] >= pezzo:
                    bins[i] -= pezzo
                    inserito = True
                    break
            if not inserito:
                bins.append(lung_comm_cm - pezzo) 
        
        sfrido_cm = sum(bins)
        riepilogo.append({
            'Diametro': f"Φ {phi}",
            'Barre 12m Intere usate': barre_intere,
            'Barre 12m da Tagliare': len(bins),
            'TOTALE BARRE DA ORDINARE': barre_intere + len(bins),
            'Sfrido Totale (metri)': round(sfrido_cm / 100, 2)
        })
        
    return pd.DataFrame(riepilogo)

def draw_arrow(msp, x, y, size, layer):
    msp.add_solid([(x, y), (x - size, y + size*0.4), (x - size, y - size*0.4)], dxfattribs={'layer': layer, 'color': 2})


# ==========================================
# 2. GENERAZIONE DEL DXF E GRAFICA
# ==========================================

def crea_armatura_definitiva(filename_dxf, file_excel, file_output_tabella, campate_mm, lung_barra_mm=12000, 
                             min_barra_mm=3000, min_sov_mm=1000, moltiplicatore_sov=50, scala_stampa=100):
    
    print("Lettura Excel e Calcolo Geometrie in corso...")
    conci_estratti, lung_totale_excel, disegna_inf, _ = carica_input_excel(file_excel)
    lung_assi = sum(campate_mm)

    barre_sup = calcola_posizionamento_barre(campate_mm, conci_estratti, 'sup', lung_barra_mm, min_barra_mm, min_sov_mm, moltiplicatore_sov, shift_sfalsamento_mm=0)
    
    # Calcolo armatura inferiore SOLO se richiesto dall'Excel
    if disegna_inf:
        sfalsamento_inferiore = lung_barra_mm / 4 
        barre_inf = calcola_posizionamento_barre(campate_mm, conci_estratti, 'inf', lung_barra_mm, min_barra_mm, min_sov_mm, moltiplicatore_sov, shift_sfalsamento_mm=sfalsamento_inferiore)

    doc = ezdxf.new('R2010', setup=True)
    msp = doc.modelspace()
    
    if 'DASHED' not in doc.linetypes:
        doc.linetypes.add('DASHED', pattern="A, 12.7, -6.35", description="Dashed __ __ __")

    doc.layers.add('01_ASSI', color=4, linetype='CENTER')
    doc.layers.add('02_CONCI', color=8)
    doc.layers.add('03_ARM_SUP', color=7)
    if disegna_inf:
        doc.layers.add('03_ARM_INF', color=3)
        doc.layers.add('04_SOVRAP_INF', color=6)
    doc.layers.add('04_SOVRAP_SUP', color=6)
    doc.layers.add('05_QUOTE_TESTI', color=2)
    doc.layers.add('06_APPOGGI', color=7) 
    doc.layers.add('07_FASI_GETTO', color=2)

    h_testo, h_frecce, gap_testo = 2.50 * scala_stampa, 2.00 * scala_stampa, 1.00 * scala_stampa
    dimstyle = doc.dimstyles.new('QUOTE_CM')
    dimstyle.dxf.dimtxt = h_testo; dimstyle.dxf.dimasz = h_frecce; dimstyle.dxf.dimgap = gap_testo; dimstyle.dxf.dimlfac = 0.1

    # Asse longitudinale
    msp.add_line((0, 0), (lung_assi, 0), dxfattribs={'layer': '01_ASSI'})
    
    # === CONCI ===
    h_conci = 15 * scala_stampa 
    for c in conci_estratti:
        msp.add_line((c['x1'], -h_conci), (c['x1'], h_conci), dxfattribs={'layer': '02_CONCI'})
        msp.add_text(c['nome'], height=h_testo*0.8, dxfattribs={'layer': '02_CONCI'}).set_placement(((c['x1']+c['x2'])/2, 5*scala_stampa), align=TextEntityAlignment.MIDDLE_CENTER)
    msp.add_line((conci_estratti[-1]['x2'], -h_conci), (conci_estratti[-1]['x2'], h_conci), dxfattribs={'layer': '02_CONCI'})

    # === CALCOLO SUPPORTO E FASI ===
    supporti_x = [0]
    x_corrente = 0
    for L in campate_mm:
        x_corrente += L
        supporti_x.append(x_corrente)

    mid_spans = [(supporti_x[i] + supporti_x[i+1])/2 for i in range(len(supporti_x)-1)]
    phases = [{'start': supporti_x[0], 'end': mid_spans[0]}]
    for i in range(len(mid_spans)-1):
        phases.append({'start': mid_spans[i], 'end': mid_spans[i+1]})
    phases.append({'start': mid_spans[-1], 'end': supporti_x[-1]})

    y_fasi_base = -120 * scala_stampa
    fasi_step = 25 * scala_stampa
    y_fasi_bottom = y_fasi_base - len(phases) * fasi_step - 10 * scala_stampa

    # === DISEGNO PILE, SPALLE E APPOGGI ===
    for idx, sx in enumerate(supporti_x):
        nome_appoggio = "SPA" if idx == 0 else "SPB" if idx == len(supporti_x) - 1 else f"P{idx}"

        th = h_testo * 1.5
        tw = len(nome_appoggio) * th * 0.9 
        pad = 2 * scala_stampa 
        y_center_box = 70 * scala_stampa
        rx1, ry1 = sx - tw/2 - pad, y_center_box - th/2 - pad
        rx2, ry2 = sx + tw/2 + pad, y_center_box + th/2 + pad
        
        msp.add_line((sx, y_fasi_bottom), (sx, ry1), dxfattribs={'layer': '06_APPOGGI', 'color': 5})

        dim_triang = 4 * scala_stampa
        y_apex = 0 
        
        msp.add_solid([(sx, y_apex), (sx - dim_triang, y_apex - dim_triang * 1.5), (sx + dim_triang, y_apex - dim_triang * 1.5)], dxfattribs={'layer': '06_APPOGGI', 'color': 7})
        msp.add_lwpolyline([(sx, y_apex), (sx - dim_triang, y_apex - dim_triang * 1.5), (sx + dim_triang, y_apex - dim_triang * 1.5), (sx, y_apex)], dxfattribs={'layer': '06_APPOGGI', 'color': 1})
        msp.add_line((sx - dim_triang * 1.5, y_apex - dim_triang * 1.5), (sx + dim_triang * 1.5, y_apex - dim_triang * 1.5), dxfattribs={'layer': '06_APPOGGI', 'color': 1})

        msp.add_lwpolyline([(rx1, ry1), (rx2, ry1), (rx2, ry2), (rx1, ry2), (rx1, ry1)], dxfattribs={'layer': '06_APPOGGI', 'color': 5})
        msp.add_text(nome_appoggio, height=th, dxfattribs={'layer': '06_APPOGGI', 'color': 1}).set_placement((sx, y_center_box), align=TextEntityAlignment.MIDDLE_CENTER)

    # === DISEGNO FASI DI GETTO ===
    for k, ph in enumerate(phases):
        fase_num = k + 1
        y_curr = y_fasi_base - k * fasi_step
        rect_h = 4 * scala_stampa
        
        hatch = msp.add_hatch(color=2, dxfattribs={'layer': '07_FASI_GETTO'})
        hatch.set_pattern_fill('ANSI31', scale=1.0)
        hatch.paths.add_polyline_path([(ph['start'], y_curr), (ph['end'], y_curr), (ph['end'], y_curr - rect_h), (ph['start'], y_curr - rect_h)], is_closed=True)
        msp.add_lwpolyline([(ph['start'], y_curr), (ph['end'], y_curr), (ph['end'], y_curr - rect_h), (ph['start'], y_curr - rect_h), (ph['start'], y_curr)], dxfattribs={'layer': '07_FASI_GETTO'})
        
        msp.add_text(f"FASE GETTO {fase_num}", height=h_testo, dxfattribs={'layer': '07_FASI_GETTO', 'color': 7}).set_placement((ph['start'] - 5*scala_stampa, y_curr - rect_h/2), align=TextEntityAlignment.MIDDLE_RIGHT)
        
        num_arrows = int((ph['end'] - ph['start']) / 4000)
        for j in range(num_arrows):
            ax = ph['start'] + 2000 + j * 4000
            draw_arrow(msp, ax, y_curr + 2*scala_stampa, 2.5*scala_stampa, '07_FASI_GETTO')
            
        if k < len(phases) - 1:
            msp.add_line((ph['end'], 0), (ph['end'], y_curr - rect_h), dxfattribs={'layer': '07_FASI_GETTO', 'linetype': 'DASHED'})

    # === FUNZIONE DI DISEGNO ARMATURE ===
    dati_tabella = []
    contatore_pos = 1

    def disegna_strato(barre, nome_strato, y_base, offset_quota_y, offset_testo_y, h_piega, layer_linea, layer_sov):
        nonlocal contatore_pos
        offset_alternanza = 3.0 * scala_stampa
        
        for i, b in enumerate(barre):
            b['y'] = y_base + (offset_alternanza if i % 2 != 0 else 0)

        for count, b in enumerate(barre):
            x1, x2, y_barra = b['x1'], b['x2'], b['y']
            lung_effettiva_mm = x2 - x1
            
            lung_totale_barra = lung_effettiva_mm + abs(h_piega) if b['type'] in ['abutment_left', 'abutment_right'] else lung_effettiva_mm
            lung_effettiva_cm = int(lung_totale_barra / 10)

            passo_mm = b['passo_cm'] * 10
            n_barre = math.ceil(b['B_mm'] / passo_mm) + 1 
            peso_unitario = (lung_totale_barra / 1000) * ((b['phi']**2) * 0.006165)
            
            dati_tabella.append({
                'Posizione': contatore_pos, 'Strato': nome_strato, 'Diametro (mm)': b['phi'], 'Passo (cm)': b['passo_cm'],
                'L Rif. B (mm)': b['B_mm'], 'Lunghezza Barra (cm)': lung_effettiva_cm, 'N. Barre': n_barre,
                'Peso Unitario (kg)': round(peso_unitario, 2), 'Peso Totale (kg)': round(peso_unitario * n_barre, 2)
            })

            sov_left = min(x2, barre[count-1]['x2']) if count > 0 and x1 < barre[count-1]['x2'] else None
            sov_right = max(x1, barre[count+1]['x1']) if count < len(barre)-1 and x2 > barre[count+1]['x1'] else None
            x_start_non = sov_left if (sov_left and sov_left > x1) else x1
            x_end_non = sov_right if (sov_right and sov_right < x2) else x2

            if sov_left and sov_left > x1: msp.add_line((x1, y_barra), (sov_left, y_barra), dxfattribs={'lineweight': 30, 'layer': layer_sov})
            if sov_right and sov_right < x2: msp.add_line((sov_right, y_barra), (x2, y_barra), dxfattribs={'lineweight': 30, 'layer': layer_sov})
            if x_end_non > x_start_non: msp.add_line((x_start_non, y_barra), (x_end_non, y_barra), dxfattribs={'lineweight': 30, 'layer': layer_linea})

            if b['type'] == 'abutment_left': msp.add_line((x1, y_barra), (x1, y_barra + h_piega), dxfattribs={'lineweight': 30, 'layer': layer_linea})
            elif b['type'] == 'abutment_right': msp.add_line((x2, y_barra), (x2, y_barra + h_piega), dxfattribs={'lineweight': 30, 'layer': layer_linea})

            msp.add_linear_dim(base=(x1, y_barra + offset_quota_y), p1=(x1, y_barra), p2=(x2, y_barra), dimstyle='QUOTE_CM', override={'layer': '05_QUOTE_TESTI'}).render()

            raggio = h_testo * 0.7
            testo_armatura = f"%%c{b['phi']}/{b['passo_cm']}  L={lung_effettiva_cm}cm"
            stima = len(testo_armatura) * h_testo * 0.7
            cx = (x1 + (lung_effettiva_mm / 2)) - ((raggio * 2) + (2 * scala_stampa) + stima) / 2 + raggio
            cy_testo = y_barra + offset_testo_y
            
            msp.add_circle((cx, cy_testo), radius=raggio, dxfattribs={'layer': '05_QUOTE_TESTI'}) 
            msp.add_text(str(contatore_pos), height=h_testo*0.8, dxfattribs={'layer': '05_QUOTE_TESTI'}).set_placement((cx, cy_testo), align=TextEntityAlignment.MIDDLE_CENTER)
            msp.add_text(testo_armatura, height=h_testo, dxfattribs={'layer': '05_QUOTE_TESTI'}).set_placement((cx + raggio + 2 * scala_stampa, cy_testo), align=TextEntityAlignment.MIDDLE_LEFT)

            contatore_pos += 1

        for i in range(len(barre) - 1):
            b1, b2 = barre[i], barre[i+1]
            x_sov_start, x_sov_end = max(b1['x1'], b2['x1']), min(b1['x2'], b2['x2'])
            if (x_sov_end - x_sov_start) > 10: 
                y_quota_sov = max(b1['y'], b2['y']) + (15 * scala_stampa if y_base > 0 else -15 * scala_stampa)
                msp.add_linear_dim(base=(x_sov_start, y_quota_sov), p1=(x_sov_start, b1['y']), p2=(x_sov_end, b2['y']), dimstyle='QUOTE_CM', override={'layer': '05_QUOTE_TESTI'}).render()

    print("Disegno Strato Superiore...")
    disegna_strato(barre_sup, 'Superiore', 35*scala_stampa, 8*scala_stampa, -4*scala_stampa, -200, '03_ARM_SUP', '04_SOVRAP_SUP')
    
    # === CONTROLLO E DISEGNO STRATO INFERIORE ===
    if disegna_inf:
        print("Disegno Strato Inferiore...")
        disegna_strato(barre_inf, 'Inferiore', -35*scala_stampa, -8*scala_stampa, 4*scala_stampa, 200, '03_ARM_INF', '04_SOVRAP_INF')
    else:
        print("Armatura Inferiore NON richiesta nell'Excel (presenza di 'None'). Salto il disegno e il calcolo.")

    doc.saveas(filename_dxf)
    print(f"\n✅ File DXF '{filename_dxf}' generato con successo!")

    # --- SALVATAGGIO EXCEL E NESTING ---
    df_distinta = pd.DataFrame(dati_tabella)
    df_nesting = esegui_nesting_sfridi(dati_tabella)

    with pd.ExcelWriter(file_output_tabella) as writer:
        totale_peso = df_distinta['Peso Totale (kg)'].sum()
        riga_totale = pd.DataFrame([{'Posizione': 'TOTALE', 'Peso Totale (kg)': round(totale_peso, 2)}])
        df_distinta_finale = pd.concat([df_distinta, riga_totale], ignore_index=True)
        df_distinta_finale.to_excel(writer, sheet_name='Distinta Ferri', index=False)
        df_nesting.to_excel(writer, sheet_name='Ottimizzazione Ordini', index=False)
        
    print(f"✅ Tabella Excel Completa salvata in '{file_output_tabella}'!")

def crea_geometria_impalcato(filename_dxf, file_excel, campate_mm, scala_stampa=100):
    print("\nGenerazione DXF Geometria con Tabella Dati in corso...")
    
    conci_estratti, _, _, _ = carica_input_excel(file_excel)
    
    doc = ezdxf.new('R2010', setup=True)
    msp = doc.modelspace()
    
    if 'DASHDOT' not in doc.linetypes:
        doc.linetypes.add('DASHDOT', pattern="A, 12.7, -6.35, 0, -6.35", description="Dash dot _ . _ . _")
        
    doc.layers.add('01_ASSI', color=4, linetype='CENTER')
    doc.layers.add('02_CONCI', color=7)
    doc.layers.add('03_QUOTE', color=2)
    doc.layers.add('04_TESTI', color=3)
    doc.layers.add('06_APPOGGI', color=7)
    
    # Sdoppiamento layer Tabella per Colori: 1=Rosso (Griglia), 2=Giallo (Testi)
    doc.layers.add('08_TAB_LINEE', color=1) 
    doc.layers.add('08_TAB_TESTI', color=2) 

    h_testo = 2.50 * scala_stampa
    dimstyle = doc.dimstyles.new('QUOTE_CM')
    dimstyle.dxf.dimtxt = h_testo; dimstyle.dxf.dimasz = 2.0 * scala_stampa; dimstyle.dxf.dimgap = 1.0 * scala_stampa; dimstyle.dxf.dimlfac = 0.1

    lung_assi = sum(campate_mm)
    msp.add_line((0, 0), (lung_assi, 0), dxfattribs={'layer': '01_ASSI'})

    y_dim_conci = 30 * scala_stampa
    y_dim_campate = 60 * scala_stampa
    y_dim_totale = 90 * scala_stampa

    # --- 1. DISEGNO CONCI E QUOTE CONCI ---
    for c in conci_estratti:
        x1, x2 = c['x1'], c['x2']
        y1_bot, y2_bot = -c['h_start'], -c['h_end']
        
        msp.add_lwpolyline([(x1, 0), (x2, 0), (x2, y2_bot), (x1, y1_bot), (x1, 0)], dxfattribs={'layer': '02_CONCI'})
        
        x_mid = (x1 + x2) / 2
        y_mid = (y1_bot + y2_bot) / 4 
        msp.add_text(c['nome'], height=h_testo*3.0, dxfattribs={'layer': '04_TESTI', 'color': 8}).set_placement((x_mid, y_mid), align=TextEntityAlignment.MIDDLE_CENTER)

        msp.add_linear_dim(base=(x1, y_dim_conci), p1=(x1, 0), p2=(x2, 0), dimstyle='QUOTE_CM', override={'layer': '03_QUOTE'}).render()
        numero_concio = c['nome'].replace('C', '')
        msp.add_text(f"CONCIO {numero_concio}", height=h_testo*0.8, dxfattribs={'layer': '04_TESTI', 'color': 2}).set_placement((x_mid, y_dim_conci - 4*scala_stampa), align=TextEntityAlignment.MIDDLE_CENTER)

    # --- 2. QUOTATURA CAMPATE E CALCOLO ASSI ---
    x_curr = 0
    supporti_x = [0]
    for L in campate_mm:
        msp.add_linear_dim(base=(x_curr, y_dim_campate), p1=(x_curr, 0), p2=(x_curr + L, 0), dimstyle='QUOTE_CM', override={'layer': '03_QUOTE'}).render()
        x_curr += L
        supporti_x.append(x_curr)

    # --- 3. QUOTATURA TOTALE ---
    msp.add_linear_dim(base=(0, y_dim_totale), p1=(0, 0), p2=(lung_assi, 0), dimstyle='QUOTE_CM', override={'layer': '03_QUOTE'}).render()

    # --- 4. DISEGNO ASSI, PILE E SPALLE ---
    y_bottom_assi = -150 * scala_stampa 
    for idx, sx in enumerate(supporti_x):
        nome_appoggio = "SPA" if idx == 0 else "SPB" if idx == len(supporti_x) - 1 else f"P{idx}"
        th = h_testo * 1.5; tw = len(nome_appoggio) * th * 0.9; pad = 2 * scala_stampa 
        y_center_box = 120 * scala_stampa
        rx1, ry1 = sx - tw/2 - pad, y_center_box - th/2 - pad
        rx2, ry2 = sx + tw/2 + pad, y_center_box + th/2 + pad
        
        # msp.add_line((sx, y_bottom_assi), (sx, ry1), dxfattribs={'layer': '06_APPOGGI', 'color': 5, 'linetype': 'DASHDOT'})
        # dim_triang = 4 * scala_stampa
        # y_apex = 0
        # msp.add_solid([(sx, y_apex), (sx - dim_triang, y_apex - dim_triang * 1.5), (sx + dim_triang, y_apex - dim_triang * 1.5)], dxfattribs={'layer': '06_APPOGGI', 'color': 7})
        # msp.add_lwpolyline([(sx, y_apex), (sx - dim_triang, y_apex - dim_triang * 1.5), (sx + dim_triang, y_apex - dim_triang * 1.5), (sx, y_apex)], dxfattribs={'layer': '06_APPOGGI', 'color': 1})
        # msp.add_line((sx - dim_triang * 1.5, y_apex - dim_triang * 1.5), (sx + dim_triang * 1.5, y_apex - dim_triang * 1.5), dxfattribs={'layer': '06_APPOGGI', 'color': 1})
        # msp.add_lwpolyline([(rx1, ry1), (rx2, ry1), (rx2, ry2), (rx1, ry2), (rx1, ry1)], dxfattribs={'layer': '06_APPOGGI', 'color': 5})
        # msp.add_text(nome_appoggio, height=th, dxfattribs={'layer': '06_APPOGGI', 'color': 1}).set_placement((sx, y_center_box), align=TextEntityAlignment.MIDDLE_CENTER)


    # ============================================================
    # --- 5. TABELLA INFERIORE (Griglia Rossa, Testi Gialli, Titoli Sospesi) ---
    # ============================================================
    
    y_tab_start = -250 * scala_stampa
    h_row = 18 * scala_stampa # Riga leggermente più alta per far respirare il testo più grande
    
    righe = [
        ("Sviluppo complessivo", "sviluppo"),
        ("Luci campate", "campate"),
        ("Concio geometrico", "nome"),
        ("Tipo concio strutturale", "strutturale"),
        ("Lunghezza concio (mm)", "L_mm"),
        ("Pioli Nelson", "pioli"),
        ("Piattabanda superiore", "psup"),
        ("Saldatura anima-pb.sup.", "sald_psup"),
        ("Altezza 1 trave (H.iniziale)", "h_start"),
        ("Altezza 2 trave (H.finale)", "h_end"),
        ("Spessore anima", "sa"),
        ("Saldatura anima-pb.inf.", "sald_pinf"),
        ("Piattabanda inferiore", "pinf")
    ]
    
    y_tab_end = y_tab_start - len(righe) * h_row
    
    # Bordo Esterno Tabella (parte da X=0, escludendo i titoli che "fluttuano")
    msp.add_lwpolyline([(0, y_tab_start), (lung_assi, y_tab_start), (lung_assi, y_tab_end), (0, y_tab_end), (0, y_tab_start)], dxfattribs={'layer': '08_TAB_LINEE'})

    # Disegno righe orizzontali e nomi intestazione fluttuanti a sinistra
    for i, r in enumerate(righe):
        y_line = y_tab_start - i * h_row
        y_text = y_line - h_row / 2
        
        # Righe di separazione orizzontali (da X=0 alla fine)
        if i > 0:
            msp.add_line((0, y_line), (lung_assi, y_line), dxfattribs={'layer': '08_TAB_LINEE'})
            
        # Nomi Righe a sinistra, fuori dalla tabella, allineati a destra (Testo Giallo Ingrandito)
        msp.add_text(r[0], height=h_testo*1.0, dxfattribs={'layer': '08_TAB_TESTI'}).set_placement((-8*scala_stampa, y_text), align=TextEntityAlignment.MIDDLE_RIGHT)

    # Disegno divisori verticali e Testi Celle
    
    # 1) Sviluppo (nessun divisore interno)
    msp.add_text(str(int(lung_assi)), height=h_testo*1.0, dxfattribs={'layer': '08_TAB_TESTI'}).set_placement((lung_assi/2, y_tab_start - h_row/2), align=TextEntityAlignment.MIDDLE_CENTER)
    
    # 2) Campate (divisori sulle pile/spalle)
    y_campate_text = y_tab_start - h_row - (h_row/2)
    x_curr = 0
    for L in campate_mm:
        msp.add_text(str(int(L)), height=h_testo*1.0, dxfattribs={'layer': '08_TAB_TESTI'}).set_placement((x_curr + L/2, y_campate_text), align=TextEntityAlignment.MIDDLE_CENTER)
        x_curr += L
        if x_curr < lung_assi: 
            msp.add_line((x_curr, y_tab_start - h_row), (x_curr, y_tab_start - 2*h_row), dxfattribs={'layer': '08_TAB_LINEE'})

    # 3+) Conci (divisori ad ogni concio dal livello 3 in giù)
    for c in conci_estratti:
        x_mid = (c['x1'] + c['x2']) / 2
        
        # Traccia la linea verticale divisoria destra del concio
        if c['x2'] < lung_assi:
            msp.add_line((c['x2'], y_tab_start - 2*h_row), (c['x2'], y_tab_end), dxfattribs={'layer': '08_TAB_LINEE'})

        # Funzione helper per piazzare il testo centrato
        def posiziona_dato(valore, indice_riga):
            y_txt = y_tab_start - (indice_riga * h_row) - (h_row / 2)
            stringa_val = str(valore).replace('.0', '') 
            if stringa_val == 'nan' or stringa_val == 'None' or stringa_val == '': stringa_val = "-"
            # Testo Giallo e più grande (h_testo*1.0)
            msp.add_text(stringa_val, height=h_testo*1.0, dxfattribs={'layer': '08_TAB_TESTI'}).set_placement((x_mid, y_txt), align=TextEntityAlignment.MIDDLE_CENTER)

        # Inserimento dati cella per cella
        posiziona_dato(c['nome'], 2)
        posiziona_dato(c['strutturale'], 3)
        posiziona_dato(c['L_mm'], 4)
        posiziona_dato(c['pioli'], 5)
        posiziona_dato(c['psup'], 6)
        posiziona_dato(c['sald_psup'], 7)
        posiziona_dato(c['h_start'], 8)
        posiziona_dato(c['h_end'], 9)
        posiziona_dato(c['sa'], 10)
        posiziona_dato(c['sald_pinf'], 11)
        posiziona_dato(c['pinf'], 12)

    doc.saveas(filename_dxf)
    print(f"✅ File Geometria e Tabella '{filename_dxf}' generati con successo!")
    
# # ==========================================
# # ESECUZIONE GLOBALE
# # ==========================================
# percorso_excel = r"C:\Users\d.gaudioso\OneDrive - Matildi+Partners\02_script\armatureLongSoletta\input_conci.xlsx"
# out_armature = r"C:\Users\d.gaudioso\OneDrive - Matildi+Partners\02_script\armatureLongSoletta\ponte_armature.dxf"
# out_geometria = r"C:\Users\d.gaudioso\OneDrive - Matildi+Partners\02_script\armatureLongSoletta\ponte_geometria.dxf"
# out_tabella = r"C:\Users\d.gaudioso\OneDrive - Matildi+Partners\02_script\armatureLongSoletta\distinta_e_ordini.xlsx"

# print("Avvio elaborazione globale...")

# # 1. LETTURA UNICA DI TUTTI I DATI DALL'EXCEL
# conci_globali, lung_totale, disegna_inf, campate_lette = carica_input_excel(percorso_excel)

# if not campate_lette:
#     print("ERRORE: Impossibile procedere senza le campate. Controlla il file Excel.")
# else:
#     # 2. ESECUZIONE MODULO ARMATURE
#     crea_armatura_definitiva(out_armature, percorso_excel, out_tabella, campate_lette)

#     # 3. ESECUZIONE MODULO GEOMETRIA
#     crea_geometria_impalcato(out_geometria, percorso_excel, campate_lette)

#     print("\n🎉 Operazione completata! Trovi i file DXF nella cartella di destinazione.")