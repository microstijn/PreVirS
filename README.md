## previr plans wur

| phase                                     | key tasks                                                                                                      |
|-------------------------------------------|----------------------------------------------------------------------------------------------------------------|
| phase 1: establish the base model (Loire) | download and run the existing Loire estuary model from the Grasso & Caillaud paper.                            |
|                                           | confirm the model reproduces the validated simulation of the estuary's physics (water, sediment, temperature). |
| phase 2: add the water quality component  | couple DELWAQ to read the output files from the hydrodynamic model.                                            |
|                                           | define viruses as a substance in DELWAQ and implement decay/transport processes using WP3 lab data.            |
| phase 3: build the shellfish module       | develop the custom module for shellfish bioaccumulation as a process within DELWAQ.                            |
|                                           | validate and calibrate the module using the oyster data from WP2.                                              |
| phase 4: test the fully integrated system | run the complete, coupled model (hydrodynamics -> DELWAQ -> shellfish).                                        |
|                                           | compare model predictions against the full dataset from WP1, WP2, and WP3.                                     |
| phase 5: run predictive scenarios         | use the validated model to simulate key events, like floods or droughts.                                       |
|                                           | generate final outputs: predictive maps and risk assessments for the Loire.                                    |
| phase 6: apply methodology to Irish site  | begin model setup for the Irish estuary (new grid, bathymetry, boundary conditions).                           |

## previr plans wur as a flowchart

```mermaid

%%{init: {'theme': 'dark', 'themeVariables': {'background': '#002060'}}}%%
graph TD
    subgraph "Partner contributions (Data & Experiments)"
        A["WP1: Field data & sampling<br>Collect water, WWTP, &<br>shellfish samples."]
        B["WP2 & WP3: Lab analysis<br>Virus quantification (dPCR),<br>diversity, and decay rates."]
    end

    subgraph "My contribution (Modeling)"
        C["1. Replicate hydrodynamic model<br>(Grasso & Caillaud, 2023).<br>This provides the physical transport<br>(currents, sediment, temp)."]
        D["2. Couple water quality model (DELWAQ)<br>Use outputs from the hydro model<br>to simulate virus transport & decay."]
        E["3. Integrate shellfish uptake module<br>Add bioaccumulation as a process<br>within the DELWAQ simulation."]
    end

    subgraph "Integrated project outcome"
        F["Calibration & validation<br>Compare the full coupled model against<br>all partner data (WP1, 2, 3)."]
        G["Scenario simulation<br>Run 'what if' scenarios (e.g., floods)<br>to generate predictive maps."]
    end

    A --> C
    A --> B
    B --> D
    B --> E
    B --> F
    C --> D
    D --> E
    E --> F
    F --> G

    style C fill:#004c99,stroke:#ccc,stroke-width:2px,color:#fff
    style D fill:#004c99,stroke:#ccc,stroke-width:2px,color:#fff
    style E fill:#004c99,stroke:#ccc,stroke-width:2px,color:#fff



```
