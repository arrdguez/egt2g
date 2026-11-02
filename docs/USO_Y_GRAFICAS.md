# Guía de ejecución y lectura de resultados

## 1. Preparar Python (Windows, Git Bash)

En la carpeta raíz del proyecto (donde está `requirements.txt`):

```bash
python --version                    # se requiere Python 3.10 o posterior
python -m venv .venv
source .venv/Scripts/activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
```

En **PowerShell**, la activación es `.venv\Scripts\Activate.ps1`; si la política de ejecución la impide, se puede invocar `.venv\Scripts\python.exe -m pip install -r requirements.txt` y después `.venv\Scripts\python.exe -m sco0 ...`, sin activar el entorno. En JetBrains se puede seleccionar `.venv\Scripts\python.exe` como intérprete del proyecto para quitar avisos de imports no resueltos.

## 2. Correr las simulaciones

```bash
python -m sco0 --config examples/one_fe.json --out results/one_fe
python -m sco0 --config examples/six_ligands.json --out results/six_ligands
python -m sco0 --config examples/two_fe.json --out results/two_fe
```

La primera es la comprobación radial simplificada, la segunda tiene **un Fe y seis N individuales** que se desplazan únicamente por su eje, y la tercera tiene **dos Fe** con acoplamiento estructural simplificado. No hace falta ejecutar las tres: cada comando escribe en la carpeta indicada por `--out`. Repetir un comando reemplaza sus resultados. `attempts.csv` puede ocupar decenas de MB y el cálculo de dos Fe puede tardar varios minutos. Los parámetros se modifican en los JSON; guarda una copia con otro nombre si quieres preservar una comparación.

## 3. Ver las gráficas, ya generadas en PNG

Abre el archivo PNG con el explorador de archivos o desde el proyecto en JetBrains. En la carpeta de cada ejecución aparecen:

| Archivo | Qué representa |
| --- | --- |
| `trajectories.png` | **Empieza aquí**: distancia Fe–N frente a `R_red`, estado de cada Fe, energía total y fracción HS. Superpone equilibrio, compresión y expansión, barridos rápidos y lentos. |
| `surfaces.png` | Energías HS y LS frente a una distancia común, diferencia `E_LS-E_HS` y separación efectiva `Δ_o(r)`. El cruce marcado supone **seis distancias iguales**. |
| `ligand_distances.png` | Solo en `six_ligands`: evolución separada de N1–N6 durante el barrido lento. N1–N4 son ecuatoriales; N5–N6 axiales. |
| `nominal_distance.png` | Fe–N frente a temperatura **nominal**, que es solo una etiqueta del protocolo geométrico. |
| `convergence.png` | Energía durante los barridos MC en el **primer** paso de una réplica; todos los pasos se guardan en `convergence.csv`. |

Por ejemplo, tras ejecutar los comandos: [gráfica principal de un Fe](../results/one_fe/trajectories.png), [seis enlaces](../results/six_ligands/ligand_distances.png) y [dos Fe](../results/two_fe/trajectories.png). Estos enlaces funcionan **localmente tras generar los resultados**; `results/` se excluye deliberadamente de Git por su tamaño. En GitHub hay que clonar y ejecutar antes de que existan las imágenes.

Si deseas mostrar una figura con Matplotlib en una ventana en vez de abrir el PNG, desde la raíz y con un entorno gráfico disponible:

```bash
python -c "import matplotlib.pyplot as plt, matplotlib.image as mpimg; plt.imshow(mpimg.imread('results/one_fe/trajectories.png')); plt.axis('off'); plt.show()"
```

## 4. Datos para análisis propio

Todos los CSV pueden abrirse en Excel, LibreOffice, JetBrains o pandas. La **fila es un paso** de una réplica en `trajectories.csv`; en `attempts.csv` **cada fila es una propuesta** (aceptada o rechazada). Los nombres con `_A` son distancias en Å, con `_eV` energías en eV. `spin1_HS=1` significa HS; `0` significa LS. `hs_fraction=0.5` significa un Fe HS y otro LS cuando hay dos centros (o una media entre réplicas en `averages.csv`).

Para ver los pasos en que un Fe cambia de estado en el barrido lento de una réplica:

```bash
python - <<'PY'
import pandas as pd
data = pd.read_csv('results/two_fe/trajectories.csv')
serie = data[(data.speed == 'slow') & (data.direction == 'cooling') & (data.replica == 0)]
print(serie.loc[serie.hs_fraction.ne(serie.hs_fraction.shift()),
                ['R_red_A', 'nominal_T_K', 'r1_A', 'r2_A', 'spin1_HS', 'spin2_HS']].to_string(index=False))
PY
```

Para analizar **únicamente las propuestas geométricas aceptadas que cambiaron el estado**:

```bash
python - <<'PY'
import pandas as pd
data = pd.read_csv('results/one_fe/attempts.csv')
cambios = data[(data.accepted == 1) & (data.spin_before_HS != data.spin_after_HS)]
print(cambios[['R_red_A', 'r_before_A', 'r_proposed_A',
               'energy_before_eV', 'energy_after_eV', 'spin_after_HS']].to_string(index=False))
PY
```

`equilibrium.csv` da la configuración de menor energía **encontrada por minimización numérica** para cada `R_red`. `spin_surfaces.csv` contiene las dos ramas (o las cuatro combinaciones para dos Fe). `averages.csv` promedia las réplicas; `config_used.json` deja constancia de los valores usados. Para seis N, `fe1_n1_A` … `fe1_n6_A` son los enlaces y `fe1_n1_x_A` … sus coordenadas relativas al Fe fijo en el origen. **La temperatura nominal no es una temperatura termodinámica ni los barridos MC representan segundos o K/min.**

## 5. Subir el código a GitHub

El proyecto debe tener ya su historial local en la rama `main`. Los datos generados, PDF y notas originales no se incluyen en el repositorio. Para crear un repositorio **privado** con GitHub CLI, primero inicia sesión interactivamente con tu propia cuenta y después, desde la raíz:

```bash
gh auth login
gh repo create egt2g-sco0 --private --source=. --remote=origin --push
```

Si ya creaste un repositorio vacío en GitHub, en lugar de `gh repo create` configura su URL real:

```bash
git remote add origin https://github.com/TU_USUARIO/TU_REPOSITORIO.git
git push -u origin main
```

Para cambios posteriores: `git add . && git commit -m "Descripción del cambio" && git push`. **No pegues tokens ni contraseñas en archivos del proyecto.**
