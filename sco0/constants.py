"""Constantes físicas y valores de referencia del campo de ligandos para Fe(II).

Unidades canónicas del proyecto:
  - energía:      cm⁻¹
  - distancia:    Å
  - temperatura:  K
  - k_B:          KB_CM1 (cm⁻¹/K)

Estas unidades son las mismas que usaba el código original `docs/prueba/main.cpp`
(energías en cm⁻¹, factor 0.695 = k_B en cm⁻¹/K).
"""

from __future__ import annotations

# Constante de Boltzmann. 1 K = 0.695038 cm⁻¹.
KB_CM1 = 0.695038

# k_B en eV/K, por si en el futuro se migra el modelo a eV.
KB_EV = 8.617333262e-5

# Parámetros de Racah del Fe(II) (repulsión electrónica e⁻–e⁻).
# En el ión libre B ≈ 1058 cm⁻¹; dentro de un complejo se reduce a ≈ 917 cm⁻¹
# (efecto nefelauxético). C ≈ 4·B. Estos valores sirven de referencia para
# conectar la linealización con el diagrama de Tanabe–Sugano completo.
RACAH_B_FREE = 1058.0
RACAH_B_COMPLEX = 917.0
RACAH_C_COMPLEX = 4.0 * RACAH_B_COMPLEX
