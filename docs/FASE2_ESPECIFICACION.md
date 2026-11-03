# Fase 2 — Especificación para implementar (modelo cooperativo → histéresis)

> **Rol:** este documento es la "inteligencia". Quien lo implemente **no debe
> razonar física ni inventar nada**: solo transcribir estas ecuaciones, estos
> parámetros y estos tests tal cual. La validación física la hará después el
> modelo `deepseek-v4-pro`.

## 1. Qué agrega esta fase

La Fase 1 (`sco0/`, ya implementada) modela **un solo centro**: el campo de
ligandos `Dq(r)` provoca el cruce HS↔LS a temperatura finita. Un centro aislado
da un cruce **suave**, sin histéresis.

La Fase 2 agrega **cooperatividad** (N centros acoplados) para obtener una
**transición de primer orden con lazo de histéresis** al calentar y enfriar.
Usa el modelo **Ising-like (campo medio + Monte Carlo)**, que es el estándar del
campo (Wajnflasz–Pick / Bousseksou) y el mismo espíritu del `project_sco` original.

**Conexión con la Fase 1:** el gap HS−LS es `Δ = 20·Dq(r) − P` (de `model.py`),
y `J` es el acoplamiento cooperativo mediado por la red elástica.

## 2. El modelo (ecuaciones exactas)

Un centro tiene espín `σ ∈ {0, 1}` (`0 = LS`, `1 = HS`). Parámetros:
`Δ` (gap, cm⁻¹), `J` (acoplamiento, cm⁻¹), `g` (degeneración HS/LS), `k_B = 0.695038 cm⁻¹/K`.

### 2.1 Ecuación de campo medio (referencia analítica)

La fracción de HS `n ∈ [0,1]` es el punto fijo de:

```
n = g / ( g + exp( (Δ − 2·J·n) / (k_B · T) ) )          (M1)
```

Para obtener el **lazo de histéresis**, resolver (M1) por iteración de punto fijo:
- **Enfriamiento:** partir de `n = g/(g+1)` a `T_max` y descender `T` reusando la `n` previa.
- **Calentamiento:** partir de `n = 0.0` a `T_min` y ascender `T` reusando la `n` previa.

### 2.2 Monte Carlo de N centros (exacto en campo medio)

Energía en función del número total de HS, `S = Σ_i σ_i`:

```
E(S) = (Δ − k_B·T·ln g)·S − (J/N)·( S·(S−1)/2 )          (M2)
```

`S·(S−1)/2` es `Σ_{i<j} σ_i σ_j`; se define `= 0` si `S < 1`. **Ojo:** (M2) depende
explícitamente de `T` (el término `−k_B·T·ln g` pliega la entropía de degeneración).
No lo "corrijas": es correcto así.

Metropolis (un solo paso = una propuesta de voltear un sitio aleatorio):
1. Elegir sitio `i` al azar.
2. `S_nuevo = S − 1` si `σ_i = 1`, o `S + 1` si `σ_i = 0`.
3. `ΔE = E(S_nuevo) − E(S)`.
4. Aceptar el volteo si `ΔE ≤ 0` o si `rand < exp(−ΔE/(k_B·T))`.
5. Actualizar `S` (y el arreglo `sigma`).

Medir `n_HS = S/N` (promedio sobre pasos de medición). Para el lazo: mismo
protocolo de enfriamiento/calentamiento que en 2.1.

## 3. Parámetros verificados (usar estos por defecto)

Valores **ya validados numéricamente por el modelo pro** para que den histéresis:

| Parámetro | Valor | Significado |
|---|---|---|
| `Delta` | 350.0 | gap HS−LS (cm⁻¹) |
| `J` | 200.0 | acoplamiento cooperativo (cm⁻¹) |
| `g` | 5.0 | degeneración HS/LS (S=2 → 2S+1=5) |
| `T_min` | 1.0 | K |
| `T_max` | 400.0 | K |
| `dT` | 1.0 | paso de temperatura (K) |
| `N` | 1000 | nº de centros (MC) |
| `mc_equil` | 50000 | pasos de termalización por temperatura |
| `mc_steps` | 50000 | pasos de medición por temperatura |
| `seed` | 2026 | semilla |

**Resultado esperado (verificado):** con estos valores el lazo va de
`T_cool ≈ 131 K` (cruce del 0.5 al enfriar) a `T_heat ≈ 266 K` (al calentar),
ancho ≈ 135 K. Condición de bistabilidad a T=0: `2·J > Δ` (aquí `400 > 350` ✓).

## 4. Archivos a crear (aditivo: NO tocar los de Fase 1)

- `sco0/cooperative.py` — todo el modelo (config, campo medio, MC, barrido, CLI).
- `examples/cooperative.json` — parámetros de la sección 3.
- `tests/test_cooperative.py` — tests de la sección 5.

Ejecutar: `python -m sco0.cooperative --config examples/cooperative.json --out results/cooperative`

## 5. Estructura de `sco0/cooperative.py` (firmas exactas)

```python
"""Fase 2: modelo cooperativo Ising-like de SCO → histéresis."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .constants import KB_CM1   # 0.695038


@dataclass(frozen=True)
class CooperativeConfig:
    Delta: float = 350.0
    J: float = 200.0
    g: float = 5.0
    T_min: float = 1.0
    T_max: float = 400.0
    dT: float = 1.0
    N: int = 1000
    mc_equil: int = 50000
    mc_steps: int = 50000
    seed: int = 2026

    def __post_init__(self): ...   # validar positivos y T_min < T_max


def mean_field(Delta: float, J: float, g: float, T: float, n0: float,
               iters: int = 2000) -> float:
    """Itera (M1) partiendo de n0; devuelve n."""
    n = n0
    for _ in range(iters):
        n = g / (g + np.exp((Delta - 2 * J * n) / (KB_CM1 * T)))
    return n


def mean_field_scan(cfg: CooperativeConfig) -> dict:
    """Devuelve {'T': array, 'n_cool': array, 'n_heat': array}."""
    # enfriar: n0 = g/(g+1), T de T_max a T_min en pasos dT
    # calentar: n0 = 0.0, T de T_min a T_max


def energy(S: int, N: int, Delta: float, J: float, g: float, T: float) -> float:
    """(M2). Σ_{i<j} σ_iσ_j = S(S-1)/2 (=0 si S<1)."""
    pairs = S * (S - 1) / 2.0 if S >= 1 else 0.0
    return (Delta - KB_CM1 * T * np.log(g)) * S - (J / N) * pairs


def mc_scan(cfg: CooperativeConfig) -> dict:
    """Metropolis sobre (M2) en ciclo térmico; devuelve n_cool y n_heat."""
    # por cada T (bajando y subiendo): termalizar mc_equil, medir mc_steps,
    # promediar n_HS = S/N. Enfriar empieza todo HS (S=N); calentar todo LS (S=0).


def main() -> None:
    """CLI: --config / --out; escribe hysteresis.csv, config_used.json y hysteresis.png."""
```

## 6. Salidas

- `hysteresis.csv` — columnas `T, n_cool_mf, n_heat_mf, n_cool_mc, n_heat_mc`.
- `hysteresis.png` — `n_HS` vs `T`; campo medio (línea) y MC (puntos), enfriamiento
  (sólido) y calentamiento (discontinuo), mismo gráfico.
- `config_used.json` — parámetros efectivos.

## 7. Tests (criterios de aceptación, valores esperados)

Escribir en `tests/test_cooperative.py` con `unittest` (sin pandas):

1. `test_high_T_limit` — `mean_field(350,200,5,10000.0, 0.5) ≈ 5/6` (±0.01).
2. `test_low_T_ls_stable` — `mean_field(350,200,5, 0.1, 0.0) ≈ 0` (< 0.01).
3. `test_energy_flip_consistency` — `energy(S±1,...) − energy(S,...)` coincide con la
   variación explícita de (M2) para S ∈ {0, 500, 999, 1000}.
4. `test_hysteresis_width` — para los defaults, `|T_cool_cruce − T_heat_cruce| > 50 K`
   (cruce = primer T donde `n` cruza 0.5).
5. `test_crossing_temperatures` — `T_cool_cruce ≈ 131 K` y `T_heat_cruce ≈ 266 K`
   (tolerancia ± 15 K).
6. `test_mc_matches_mean_field` — con `N=1000`, el MC y el campo medio coinciden en
   `n_HS` dentro de ±0.05 en un par de temperaturas lejos del salto (p. ej. 100 K y 320 K).
7. `test_reproducible` — `mc_scan` con la misma semilla da el mismo resultado.

## 8. Reglas duras (no violar)

- No cambiar las ecuaciones (M1)/(M2) ni los parámetros de la sección 3.
- No tocar los archivos de Fase 1 (`config.py`, `model.py`, `metropolis.py`, etc.).
- Unidades fijas: energía en cm⁻¹, temperatura en K, `k_B = 0.695038`.
- No añadir física nueva (presión, red explícita, movimientos de distancia). Eso es Fase 3+.
