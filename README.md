# OpenFusionSim
Open-Source Nuclear Fusion reactor simulator game

## Scope

OpenFusionSim is a nuclear fusion reactor simulator game that allows users to control a virtual fusion reactor. The game will simulate the physics of nuclear fusion and provide an interactive experience for players to learn about the principles of fusion energy.

## Features
- Realistic physics simulation of nuclear fusion reactions
- Interactive control panel for managing reactor parameters
- Educational content about fusion energy and its applications
- Multiple levels of difficulty to cater to different skill levels

## Plasma physics model

The 0-D burning-plasma core includes:

- **Two-temperature plasma** — separate electron/ion energy channels coupled by
  Coulomb collisional equilibration; alpha heating split e/i via the Stix
  slowing-down integral; ECRH/LHCD heat electrons, ICRH heats ions, NBI splits
  50/50; radiation cools the electrons only.
- **L/H confinement regimes** — Martin-2008 power threshold with hysteresis;
  H98 = 0.55 in L-mode vs 1.0 in H-mode. You must pay the L→H power toll
  before burn-grade confinement becomes available.
- **Physical ohmic heating** — P = I²R with neoclassical (trapped-particle
  corrected) Spitzer resistivity, and a published loop voltage.
- **Impurity line radiation** — coronal cooling curves (low-Z wall material
  plus a tungsten trace). A cold, dirty plasma faces a genuine radiation
  barrier at start-up: fuel before burn-through and it becomes insurmountable.
- **Sawtooth oscillations** — q₀ < 1 triggers Kadomtsev-style crashes that
  flatten the core current profile and expel a slice of stored energy.
- **Uckan q95** — elongation, triangularity and aspect-ratio corrections, so
  the plasma shape knobs genuinely move the edge safety factor.
- **Persistence-based Greenwald limit** — density-limit disruptions develop
  over seconds (radiative collapse), not as an instant trip wire.
- **Type-I ELMs** — the price of H-mode: the edge pedestal periodically
  crashes, firing MJ-scale heat pulses at the divertor tiles. Pace them by
  running the pellet injector above the natural ELM frequency (ΔW·f ≈ const),
  or watch your strike points flash toward the tungsten melt limit.
- **CS volt-second budget** — the central solenoid has a finite transformer
  flux swing (~90 Wb). The burn consumes it at V_loop × (inductive current
  fraction); when it's gone, only non-inductive current (NBI/LHCD drive +
  bootstrap) can hold I_p. This is the pulse-length limit of a conventional
  tokamak — steady state must be earned with current drive.
- **Runaway electron beams** — an unmitigated fast current quench above
  ~4 MA avalanche-converts plasma current into a relativistic electron beam
  that strikes the first wall with hundreds of MJ. Keeping MGI/SPI armed
  (auto-MGI) suppresses the avalanche — that's why the mitigation systems
  exist.
- **NBI shine-through interlock** — the 1 MeV beams are blocked below
  n_e = 1.5×10¹⁹ m⁻³ (they'd pass through a thin plasma and melt the far
  wall). Build density before bringing on the beams.
- **Neutron wall damage** — first-wall fluence and displacement damage (dpa)
  accumulate over the campaign; the ITER-class wall is rated for ~3 dpa.
- Temperature is **earned, not dialled**: there is no T_e actuator. Heat the
  plasma through the H&CD systems, exactly like a real machine.

The **OP SPACE** tab shows a live POPCON operating-space diagram — the
auxiliary power required for steady state across the (T_i, n_e) plane, with
ignition and Q=10 contours, the Greenwald/Troyon/L-H limits, and a breadcrumb
trail of your discharge — alongside the full power ledger and
limit-proximity bars.

## Installation
To install OpenFusionSim, follow these steps:
1. Clone the repository: `git clone https://github.com/PolskaKrowa/OpenFusionSim.git`
2. Navigate to the project directory: `cd OpenFusionSim`
3. Build the project using CMake: `cmake -B build && cmake --build build -j$(nproc)`
4. Run the game: `./build/bin/fusionsim`

## Contributing
Contributions to OpenFusionSim are welcome! If you have an idea for a new feature or want to fix a bug, please submit a pull request. Make sure to follow the coding standards and include tests for any new functionality.
