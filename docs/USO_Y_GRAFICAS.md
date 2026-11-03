# Guía de ejecución y lectura de resultados

## 1. Preparar Python

En la raíz del proyecto (donde está `requirements.txt`):

```bash
python --version                    # se requiere Python 3.10 o posterior
python -m venv .venv
source .venv/bin/activate           # Linux/macOS (en Windows: .venv/Scripts/activate)
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
```

Solo se necesitan `numpy` y `matplotlib`.

## 2. Correr las simulaciones

```bash
# Barrido de distancia (el campo de ligandos provoca el cruce HS↔LS)
python -m sco0 --config examples/one_fe.json --out results/one_fe

# Barrido de temperatura a distancia fija (cruce suave, no cooperativo)
python -m sco0 --config examples/one_fe_temperature.json --out results/one_fe_temperature
```

Cada comando escribe en la carpeta indicada por `--out`; repetirlo reemplaza
los resultados. Los parámetros se editan en los JSON (ver claves en
`examples/one_fe.json` o en `sco0/config.py`).

## 3. Ver las gráficas

En cada carpeta de resultados aparecen dos PNG:

| Archivo | Qué representa |
| --- | --- |
| `surfaces.png` | `E_HS(r)`, `E_LS(r)`, `Dq(r)` y `Δo(r)=10·Dq(r)` frente a la distancia, con la distancia de cruce `r_c` marcada. Muestra cómo el campo de ligandos crece al acortarse el enlace y cruza las dos superficies. |
| `crossover.png` | **Empieza aquí**: fracción de HS frente a `r` (o frente a `T`). Superpone la solución analítica (línea) y el Monte Carlo a temperatura finita (puntos). |

Abrir las imágenes con el explorador de archivos o, en un entorno gráfico:

```bash
python -c "import matplotlib.pyplot as plt, matplotlib.image as mpimg; plt.imshow(mpimg.imread('results/one_fe/crossover.png')); plt.axis('off'); plt.show()"
```

## 4. Datos para análisis propio

Los CSV se abren en Excel/LibreOffice/pandas. Cada fila es un paso del barrido:

- `equilibrium.csv`: referencia analítica por paso — `r_A` (Å), `T_K`, `Dq_cm1`,
  `delta_o_cm1`, `gap_cm1` (= `E_HS−E_LS`) y `hs_fraction` (fracción HS exacta).
- `monte_carlo.csv`: fracción HS muestreada por Metropolis (promedio entre réplicas).
- `config_used.json`: configuración efectiva usada.

Las energías van en **cm⁻¹** y las distancias en **Å**. `hs_fraction=0` es todo
LS y `hs_fraction=1` es todo HS.

## 5. Subir el código a GitHub

```bash
git add .
git commit -m "Descripción del cambio"
git push
```

No pegues tokens ni contraseñas en archivos del proyecto.
