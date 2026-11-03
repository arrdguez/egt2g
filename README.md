# SCO — transición de espín Fe(II) por campo de ligandos

Modelo mínimo de la transición de espín (spin crossover) de Fe(II). Captura el
**mecanismo de raíz** del fenómeno: el **campo de ligandos `Dq(r)`**, que depende
de la distancia Fe–N, es lo que provoca el cruce entre **alto espín (HS)** y
**bajo espín (LS)**.

> **Propósito de la modelación (el espíritu del código original `docs/prueba/main.cpp`):**
> variar la distancia Fe–N → cambiar el campo de ligandos `Dq` → y ver cómo
> **es el campo de ligandos el que provoca el cruce HS↔LS en su raíz**.
> Sobre esta base se añadirán luego los factores que inciden sobre `Dq`
> (temperatura, presión, cooperatividad de red, sustitución de ligandos).

## El modelo

Un solo centro Fe con distancia `r` (Å). El campo de ligandos decae con la
distancia (modelo de cargas puntuales, `10·Dq ∝ r⁻⁵`):

```
Dq(r) = Dq0 · (r0/r)⁵          [cm⁻¹]
```

Energías electrónicas (linealización del diagrama de Tanabe–Sugano para d⁶):

```
E_HS(r) = −4·Dq(r)             (t₂g⁴ e_g²)
E_LS(r) = −24·Dq(r) + P        (t₂g⁶, P = energía de apareamiento)
```

El cruce HS↔LS ocurre en `r_c = r0·(20·Dq0/P)^(1/5)`, donde `ΔE = E_HS − E_LS
= 20·Dq − P = 0`. A distancia larga (campo débil) domina **HS**; a distancia
corta (campo fuerte) domina **LS** (enlace Fe–N más corto, contracción de red).

**Unidades:** energía en cm⁻¹, distancia en Å, temperatura en K, con
`k_B = 0.695 cm⁻¹/K`. Valores de referencia (Fe²⁺): `Dq0 ≈ 1500 cm⁻¹`,
`B(Racah) ≈ 917–1050 cm⁻¹`, `P ≈ 20000 cm⁻¹`.

## Temperatura finita (Monte Carlo)

El espín **se propone y se acepta** con Boltzmann `min(1, exp(−ΔE/k_B T))`
(`sco0/metropolis.py`). A diferencia del prototipo anterior (que solo aceptaba
movimientos cuesta abajo a 0 K), esto recupera el ensemble térmico real del
código original. La fracción de HS tiene solución analítica exacta (sistema de
dos niveles), que se usa como referencia:

```
n_HS(r, T) = 1 / (1 + exp(ΔE / k_B T))
```

## Instalación y ejecución

```bash
python -m venv .venv
source .venv/bin/activate            # Linux/macOS
python -m pip install -r requirements.txt

# pruebas
python -m unittest discover -s tests -v

# barrido de distancia (campo de ligandos → cruce)
python -m sco0 --config examples/one_fe.json --out results/one_fe

# barrido de temperatura a distancia fija (cruce suave, no cooperativo)
python -m sco0 --config examples/one_fe_temperature.json --out results/one_fe_temperature
```

## Salidas (`results/<out>/`)

- `surfaces.png` — `E_HS(r)`, `E_LS(r)`, `Dq(r)` y `Δo(r)` con el cruce `r_c`.
- `crossover.png` — fracción HS vs `r` (o vs `T`): analítico frente a Monte Carlo.
- `equilibrium.csv` — referencia analítica por paso (con `Dq`, `Δo`, `gap`, `n_HS`).
- `monte_carlo.csv` — fracción HS muestreada por Metropolis.
- `config_used.json` — configuración efectiva.

## Qué NO modela (todavía)

Un centro aislado da un **cruce suave**, no una transición de primer orden con
histéresis. La histéresis exige **cooperatividad** (red elástica acoplada),
que es la siguiente fase (ver `docs/GUIA_PARA_CONTINUAR.pdf`). Ver también la
memoria del proyecto en `docs/MEMORIA_DEL_PROYECTO_ORIGINAL.pdf`.
