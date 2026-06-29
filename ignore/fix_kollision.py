# =============================================================================
# FIX für spl_simulation_v5 copy.ipynb
# =============================================================================
#
# Kopiere den Code der jeweiligen Zelle in die entsprechende Notebook-Zelle.
# Jeder Abschnitt ist klar mit "ZELLE X" markiert.
#
# Es gibt 3 Fixes:
#   1. Simulation.simulate()  — dynamisches Hinzufügen neuer Körper nach Kollision
#   2. simuliere_schuss()     — rollenbasierte Trajektorien-Aufzeichnung
#   3. animate_columbiade()   — korrektes Handling nach Kollision (max statt min)
# =============================================================================


# ─────────────────────────────────────────────────────────────────────────────
# ZELLE: Simulation-Klasse (Zelle 3 im Notebook)
# Ersetze die GESAMTE Zelle mit diesem Code.
# ─────────────────────────────────────────────────────────────────────────────

class Simulation:

    # N-Körper-Gravitations-Simulation mit konfigurierbarem Zeitschritt.
    # Erkennt und behandelt inelastische Kollisionen.

    def __init__(self, koerper, dt_s=60.0):
        
        # Parameter:

        # koerper : list[Koerper] – Liste der zu simulierenden Körper
        # dt_s    : float         – Zeitschritt in Sekunden
        
        if not koerper:
            raise ValueError("Mindestens ein Körper muss angegeben werden")
        if dt_s <= 0:
            raise ValueError(f"Zeitschritt muss positiv sein, erhalten: {dt_s}")

        self.koerper = list(koerper)
        self.dt = dt_s          # s
        self.zeit = 0.0         # s

    def _beschleunigungen(self):
        # Berechnet die Gravitationsbeschleunigung aller Körper (m/s²)
        n = len(self.koerper)
        acc = [np.zeros(3) for _ in range(n)]

        for i in range(n):
            for j in range(i + 1, n):
                r_vec = self.koerper[j].pos - self.koerper[i].pos
                r = np.linalg.norm(r_vec)
                if r == 0:
                    continue
                # Gravitationsbeschleunigung: a = G*m / r²  (Richtung: r_vec/r)
                g_faktor = G / r**2
                a_i = g_faktor * self.koerper[j].masse * (r_vec / r)
                a_j = g_faktor * self.koerper[i].masse * (-r_vec / r)
                acc[i] += a_i
                acc[j] += a_j

        return acc

    def _kollisionen_erkennen(self):
        # Prüft alle Körperpaare auf Berührung und führt ggf. inelastische Kollision durch
        i = 0
        while i < len(self.koerper):
            j = i + 1
            while j < len(self.koerper):
                k1 = self.koerper[i]
                k2 = self.koerper[j]
                abstand = np.linalg.norm(k2.pos - k1.pos)
                beruehrungsabstand = k1.radius + k2.radius

                if abstand <= beruehrungsabstand:
                    # 100% inelastische Kollision: Impulserhaltung
                    neue_masse = k1.masse + k2.masse
                    neue_vel = (k1.masse * k1.vel + k2.masse * k2.vel) / neue_masse
                    neue_pos = (k1.masse * k1.pos + k2.masse * k2.pos) / neue_masse
                    # Neuer Radius aus Massenerhaltung (ρ = const)
                    neuer_radius = (k1.radius**3 + k2.radius**3)**(1/3)

                    merged = Koerper(
                        name=f"{k1.name} + {k2.name}",
                        masse=neue_masse,
                        durchmesser_km=neuer_radius * 2 / 1e3,
                        position_km=neue_pos / 1e3,
                        geschwindigkeit_km_s=neue_vel / 1e3,
                        farbe=k1.farbe
                    )
                    #print(f"Kollision: {merged.name} [t={self.zeit/86400:.2f} d]")
                    self.koerper.pop(j)
                    self.koerper[i] = merged
                else:
                    j += 1
            i += 1

    def schritt(self):
        # Führt einen Euler-Integrationsschritt durch
        acc = self._beschleunigungen()
        for i, k in enumerate(self.koerper):
            k.vel += acc[i] * self.dt
            k.pos += k.vel * self.dt
        self.zeit += self.dt
        self._kollisionen_erkennen()

    def simulate(self, dauer_tage, schritte_pro_tag=None):

        # Simuliert für eine gegebene Dauer und speichert Trajektorien.

        # Parameter:

        # dauer_tage      : float – Simulationsdauer in Tagen
        # schritte_pro_tag: int   – Anzahl Zeitschritte pro Tag (überschreibt dt)
 
        # Rückgaben:

        # trajektorien : dict {name: np.ndarray (N,3)} – Positionen in km
        # zeiten       : np.ndarray – Zeiten in Tagen

        if schritte_pro_tag is not None:
            self.dt = 86400.0 / schritte_pro_tag

        gesamt_schritte = int(dauer_tage * 86400 / self.dt)
        trajektorien = {k.name: [] for k in self.koerper}
        zeiten = []

        for _ in range(gesamt_schritte):
            # Positionen speichern
            for k in self.koerper:
                # FIX: Neue Körpernamen nach Kollision dynamisch hinzufügen
                if k.name not in trajektorien:
                    trajektorien[k.name] = []
                trajektorien[k.name].append(k.pos.copy())
            zeiten.append(self.zeit)
            self.schritt()

        # Konvertierung zu numpy-Arrays in km
        for name in trajektorien:
            if len(trajektorien[name]) > 0:
                trajektorien[name] = np.array(trajektorien[name]) / 1e3  # km
            else:
                trajektorien[name] = np.empty((0, 3))
        zeiten = np.array(zeiten) / 86400  # Tage
        return trajektorien, zeiten


# ─────────────────────────────────────────────────────────────────────────────
# ZELLE: Hilfsfunktionen (Zelle 8 im Notebook — die mit "Hilfsfunktionen geladen")
# Ersetze die GESAMTE Zelle mit diesem Code.
# ─────────────────────────────────────────────────────────────────────────────

GESCHOSS_MASSE       = 1e4 # kg
GESCHOSS_DURCHMESSER = 0.01 # km

def winkel_mond(mond_pos_km):
    # Gibt den aktuellen Winkel des Mondes in der X/Y-Ebene zurück (Bogenmaß).
    return np.arctan2(mond_pos_km[1], mond_pos_km[0])

def erstelle_geschoss(abschuss_winkel_rad, geschwindigkeit_km_s):

    # Erzeugt ein Geschoss an der Erdoberfläche.

    # Parameter:

    # abschuss_winkel_rad : float – Winkel in der X/Y-Ebene (Bogenmaß)
    # geschwindigkeit_km_s: float – Mündungsgeschwindigkeit in km/s

    # Rückgaben:

    # Koerper – das Geschoss

    r_erde_km = ERDE_DURCHMESSER / 2   # km
    richtung = np.array([np.cos(abschuss_winkel_rad),
                         np.sin(abschuss_winkel_rad),
                         0.0])
    pos_km = richtung * r_erde_km
    vel_km_s = richtung * geschwindigkeit_km_s
    return Koerper('Geschoss', GESCHOSS_MASSE, GESCHOSS_DURCHMESSER,
                   position_km=pos_km,
                   geschwindigkeit_km_s=vel_km_s,
                   farbe='orange')

def _rolle_bestimmen(name):
    # Bestimmt die "Rolle" eines Körpers anhand seines Namens.
    # Nach einer Kollision heißt ein Körper z.B. "Mond + Geschoss".
    # Diese Funktion ordnet ihn der Rolle "Mond" zu, damit die
    # Trajektorie nahtlos weitergeführt wird.
    #
    # Priorität: Erde > Mond > Geschoss
    # (weil "Erde + Geschoss" als Erde weitergeführt werden soll,
    #  und "Mond + Geschoss" als Mond)
    if 'Erde' in name:
        return 'Erde'
    if 'Mond' in name:
        return 'Mond'
    if name == 'Geschoss':
        return 'Geschoss'
    return name

def simuliere_schuss(geschwindigkeit_km_s, vorlauf_grad=0.0, dauer_tage=10, dt_s=60):
    
    # Simuliert einen Kanonenschuss
    
    # Parameter:
    
    # geschwindigkeit_km_s : float – Mündungsgeschwindigkeit
    # vorlauf_grad         : float – Winkelvorlauf relativ zur Mondposition (Grad)
    # dauer_tage           : float – Simulationsdauer
    # dt_s                 : float – Zeitschritt in Sekunden
    
    # Rückgaben:
    
    # traj   : dict  – Trajektorien aller Körper in km
    # zeiten : array – Zeiten in Tagen
    # treffer: bool  – True wenn Geschoss den Mond trifft
    
    erde, mond = erstelle_erde_mond()

    # Abschusswinkel: Mondrichtung + Vorlaufwinkel
    phi_mond = winkel_mond(mond.pos / 1e3)
    phi_schuss = phi_mond + np.deg2rad(vorlauf_grad)
    geschoss = erstelle_geschoss(phi_schuss, geschwindigkeit_km_s)

    sim = Simulation([erde, mond, geschoss], dt_s=dt_s)
    treffer = False

    gesamt_schritte = int(dauer_tage * 86400 / dt_s)
    # FIX: traj wird mit Rollen-Keys initialisiert
    traj = {'Erde': [], 'Mond': [], 'Geschoss': []}
    zeiten = []

    for _ in range(gesamt_schritte):
        # FIX: Positionen nach "Rolle" speichern, nicht nach exaktem Namen.
        # Nach einer Kollision heißt der Körper z.B. "Mond + Geschoss",
        # wird aber weiterhin unter der Rolle "Mond" gespeichert.
        for k in sim.koerper:
            rolle = _rolle_bestimmen(k.name)
            if rolle in traj:
                traj[rolle].append(k.pos.copy())

        zeiten.append(sim.zeit)
        sim.schritt()

        # Treffer prüfen: Geschoss verschwunden (mit Mond verschmolzen)?
        namen = [k.name for k in sim.koerper]
        if 'Geschoss' not in namen:
            treffer = any('Mond' in n and 'Geschoss' in n for n in namen) or \
                      any('Geschoss' in n and 'Mond' in n for n in namen) or \
                      len(sim.koerper) < 3
            # War die Kollision mit dem Mond (nicht der Erde)?
            for k in sim.koerper:
                if 'Mond' in k.name and 'Geschoss' in k.name:
                    treffer = True
                if 'Erde' in k.name and 'Geschoss' in k.name:
                    treffer = False
            

    for name in traj:
        if len(traj[name]) > 0:
            traj[name] = np.array(traj[name]) / 1e3 # km
        else:
            traj[name] = np.empty((0, 3))
    zeiten = np.array(zeiten) / 86400 # Tage
    return traj, zeiten, treffer

print("Hilfsfunktionen geladen")


# ─────────────────────────────────────────────────────────────────────────────
# ZELLE: Animation (Zelle 10 im Notebook — animate_columbiade)
# Ersetze die GESAMTE Zelle mit diesem Code.
# ─────────────────────────────────────────────────────────────────────────────

def animate_columbiade(traj, zeiten, anim_v, anim_delta, treffer_anim, n_frames=200):

    # Animiert das Erde Mond Geschoss System.

    # FIX: max() statt min() — damit die Animation auch nach der Kollision
    # weiterläuft (Geschoss-Trajektorie ist kürzer als Mond/Erde)
    n_data = max(len(v) for v in traj.values())
    schritt = max(1, n_data // n_frames)
    idx = np.arange(0, n_data, schritt)

    fig, ax = plt.subplots(figsize=(6, 6))

    if anim_v == 9:
        lim = 30000
        VIS_R_GESCHOSS = 600
        VIS_R_ERDE    = 6370
    elif anim_v == 7:
        lim = 30000
        VIS_R_GESCHOSS = 600
        VIS_R_ERDE    = 6370
    else:
        lim = 500000
        VIS_R_GESCHOSS = 6000
        VIS_R_ERDE    = 15000

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect('equal')
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')
    ax.set_xlabel('X [km]', color='black')
    ax.set_ylabel('Y [km]', color='black')
    ax.tick_params(colors='black')
    ax.grid(True, alpha=0.1, color='black')

    spur_mond,    = ax.plot([], [], color='grey', lw=0.5, alpha=0.5, zorder=1)
    spur_geschoss,= ax.plot([], [], color='orange',  lw=1.2, alpha=0.8, zorder=2)

    VIS_R_MOND    = 4050

    erde_k = plt.Circle((0, 0), VIS_R_ERDE, color='blue', zorder=3)
    mond_k = plt.Circle((MOND_ABSTAND, 0), VIS_R_MOND, color='grey', zorder=3)
    gesc_k = plt.Circle((0, 0), VIS_R_GESCHOSS, color='orange', zorder=5)
    ax.add_patch(erde_k)
    ax.add_patch(mond_k)
    ax.add_patch(gesc_k)

    zeittext = ax.text(0.02, 0.95, '', transform=ax.transAxes, color='black', fontsize=9)
    treffertext = ax.text(0.5, 0.04, '', transform=ax.transAxes, color='lime', fontsize=12, ha='center', fontweight='bold')

    mond_traj_a = traj.get('Mond', traj.get(list(traj.keys())[1]))
    gesc_traj_a = traj.get('Geschoss', None)
    erde_traj_a = traj.get('Erde', None)

    def update(frame):
        i = idx[frame]

        # FIX: Bounds-Check für Mond-Trajektorie
        if i < len(mond_traj_a):
            mx, my = mond_traj_a[i, 0], mond_traj_a[i, 1]
            mond_k.center = (mx, my)
            spur_mond.set_data(mond_traj_a[max(0, i - 200 * schritt):i, 0], mond_traj_a[max(0, i - 200 * schritt):i, 1])

        # FIX: Bounds-Check für Erde-Trajektorie
        if erde_traj_a is not None and i < len(erde_traj_a):
            erde_k.center = (erde_traj_a[i, 0], erde_traj_a[i, 1])

        if gesc_traj_a is not None and i < len(gesc_traj_a):
            gx, gy = gesc_traj_a[i, 0], gesc_traj_a[i, 1]
            gesc_k.center = (gx, gy)
            gesc_k.set_visible(True)
            spur_geschoss.set_data(gesc_traj_a[:i, 0], gesc_traj_a[:i, 1])
        else:
            gesc_k.set_visible(False)
            if treffer_anim:
                treffertext.set_text('Treffer')

        t_tage = zeiten[i] if i < len(zeiten) else zeiten[-1]
        zeittext.set_text(f't = {t_tage:.2f} d  |  Geschwindigkeit={anim_v} km/s  Winkel={anim_delta:.0f}°')
        return erde_k, mond_k, gesc_k, spur_mond, spur_geschoss, zeittext, treffertext

    anim = FuncAnimation(fig, update, frames=len(idx), interval=35, blit=True)
    plt.close()
    return anim
