# Referencias del proyecto (SCO / transición de espín)

Material bibliográfico copiado desde `~/Dropbox/Book` y `~/Dropbox/src/prueba`
para tener todo el soporte teórico del doctorado junto al código.

> **Nota:** los PDF/DJVU están **excluidos de Git** (`*.pdf` en `.gitignore`).
> Son referencias locales de lectura; si en algún momento quieres versionarlas,
> tendrás que quitar `*.pdf` del `.gitignore` (no recomendado: son archivos grandes
> y con copyright).

---

## Libros (`libros/`)

| Archivo | Referencia | Por qué es relevante |
|---|---|---|
| `Lever_Inorganic_Electronic_Spectroscopy.pdf` | A. B. P. Lever, *Inorganic Electronic Spectroscopy* (2.ª ed., Elsevier) | **La fuente de Dq, Δo, Tanabe–Sugano y parámetros de Racah A/B/C.** De aquí salen los diagramas de niveles d y el cruce HS/LS que usaste en `prueba/main.cpp`. |
| `Kahn_Molecular_Magnetism_1993.pdf` | O. Kahn, *Molecular Magnetism*, VCH (1993) | La referencia clásica de **transición de espín (SCO)** y magnetismo molecular; base conceptual del proyecto. |
| `Theoretical_Foundations_of_Molecular_Magnetism_1999.pdf` | R. Boča, *Theoretical Foundations of Molecular Magnetism*, Elsevier (1999) | Fundamentos formales de magnetoquímica; niveles de energía, acoplamientos, ZFS. |
| `Gray_Electrons_and_Chemical_Bonding.pdf` | H. B. Gray, *Electrons and Chemical Bonding* (1965) | Teoría de campo de ligandos / OM de forma pedagógica. |
| `Kettle_Physical_Inorganic_Chemistry.pdf` | S. F. A. Kettle, *Physical Inorganic Chemistry* | Campo de ligandos y espectroscopía inorgánica. |
| `Functional_Molecular_Materials_Introductory_Textbook.pdf` | *Functional Molecular Materials: An Introductory Textbook* | Materiales moleculares funcionales (SCO, imanes moleculares). |
| `Metal_organic_and_organic_molecular_magnets.djvu` | *Metal-Organic and Organic Molecular Magnets* (RSC) | Imanes moleculares metálico-orgánicos. |
| `Molecular_Magnets_Physics_and_Applications.pdf` | J. Bartolomé, F. Luis, J. F. Fernández (eds.), *Molecular Magnets: Physics and Applications*, Springer (NanoScience and Technology) | Física y aplicaciones de imanes moleculares. |
| `Multifunctional_Molecular_Magnets_Octacyanidometalates.pdf` | K. Imoto, *Multifunctional Molecular Magnets Based on Octacyanidometalates*, Springer Theses | Análogos cianometalato (misma familia que los azules de Prusia / Fe(py)₂[Ni(CN)₄]). |
| `NMR_MRI_SR_Mossbauer_Spectroscopies_Molecular_Magnets.pdf` | *NMR–MRI, μSR and Mössbauer Spectroscopies in Molecular Magnets* | Técnicas (Mössbauer sobre todo) para caracterizar SCO. |

## Artículos (`articulos/`)

| Archivo | Referencia | Por qué es relevante |
|---|---|---|
| `Levchenko_2014_SCO_under_pressure_LowTempPhys_40_571.pdf` | G. G. Levchenko, A. V. Khristov, V. N. Varyukhin, *Spin crossover in iron(II)-containing complex compounds under a pressure*, Low Temp. Phys. **40**, 571 (2014). DOI 10.1063/1.4891445 | **El paper del que salió `prueba/main.cpp`** (SCO en función de Dq, con los valores numéricos Dq/B y el cruce). |
| `DFT_calculations_spin_crossover.pdf` | H. Paulsen, A. X. Trautwein, *Density Functional Theory Calculations for Spin Crossover Complexes*, Top. Curr. Chem. **235**, 197 (2004) | Vía DFT para calcular complejos SCO (el nivel más "microscópico" del problema). |
| `Zero_Field_Splitting_and_G_factors.pdf` | T. Bodenstein, K. Fink, *Ab initio calculations on the magnetic properties of transition metal complexes*, AIP Conf. Proc. **1702**, 090011 (2015) | Cálculo ab initio de propiedades magnéticas (D, g). |
| `JChemEd_Molecular_Integrals_2018.pdf` | K. V. Murphy, J. M. Turney, H. F. Schaefer III, *Student-Friendly Guide to Molecular Integrals*, J. Chem. Educ. **95**, 1572 (2018) | Integrales moleculares (utilidad si se avanza hacia métodos de estructura electrónica). |

---

## Material complementario (no copiado, disponible en `~/Dropbox/`)

- **Tesis**: `~/Dropbox/PhD/thesis/` (`Chapter_1.docx` … `Chapter_4.docx`, `Chapter_3_Results.docx`).
- **Datos experimentales**: `~/Dropbox/PhD/2018/Semestre-I/SCO_DATA/` (ZFCFC, MvsH, AC en piridinas sustituidas).
- **Códigos históricos**: `~/Dropbox/src/` (`DMF.beta.*`, `one_D`, `project_sco`, `prueba`, `project_one`).
- Más libros en `~/Dropbox/Book/` (Shriver & Atkins, Atkins Physical Chemistry, etc.).
