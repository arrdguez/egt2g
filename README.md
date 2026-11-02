# SCO0: spin crossover electrónico-estructural a 0 K

Prototipo **ilustrativo** de Fe(II) HS/LS con 1 o 2 centros. Usa Python 3.10+; no hace falta compilador C++. La geometría Fe–N cambia según la superficie electrónica y la distancia que intenta imponer la red. **No modela magnetismo, temperatura termodinámica, vibraciones ni cinética experimental.**

**¿Cómo lo corro y cómo veo las gráficas?** Consulta la [guía de uso y visualización](docs/USO_Y_GRAFICAS.md). Cada ejecución crea automáticamente CSV y PNG en `results/`; abre primero `results/one_fe/trajectories.png` después de ejecutar el ejemplo. Los resultados son regenerables y no se suben al repositorio Git.

## Instalación y ejecución (desde esta carpeta)

```bash
python -m venv .venv
source .venv/Scripts/activate   # Git Bash en Windows
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python -m sco0 --config examples/one_fe.json --out results/one_fe
python -m sco0 --config examples/six_ligands.json --out results/six_ligands
python -m sco0 --config examples/two_fe.json --out results/two_fe
```

Si Python ya tiene las dependencias, basta ejecutar los tres últimos comandos. Para ayuda: `python -m sco0 --help`. En PowerShell la activación del entorno virtual es distinta; también se puede llamar directamente a `.venv/Scripts/python.exe`.

## Modelo

Para cada centro, `x=r/r_ref` y `E_s(r)=E_scale*(1/(9*x**9) - 0.02*ligand_strength*rho_s/x**5 - 1/x)+P_s`. Para HS: `rho=2, P_s=0`; para LS: `rho=12, P_s=P_LS`. El potencial de Jørgensen es **ilustrativo** y originalmente no constituye una parametrización de Fe(py)₂[Ni(CN)₄]. Las energías están en eV, las distancias en Å y las rigideces en eV/Å².

En `geometry="radial"` hay una sola coordenada efectiva `r` por Fe. En `geometry="six_radial"` hay un Fe **fijo** en el origen y seis N sobre los semiejes `+x, -x, +y, -y, +z, -z`; cada posición tiene **su propio radio independiente**. N1–N4 son ecuatoriales y N5–N6 axiales. Solo **un** nitrógeno se mueve por intento, a lo largo de su eje; no hay movimiento angular ni simulación de orbitales o interacciones N–N. La energía electrónica es la **media de los seis términos** de Jørgensen más `P_LS` una sola vez. Así recupera exactamente el modelo radial cuando las seis distancias son iguales. `Δ_o = 0.1*E_scale*ligand_strength*mean((r_ref/r_j)**5)`; la distancia media *no* se inserta en la potencia quinta. La definición de Δ_o para geometrías distorsionadas es una aproximación isotrópica ilustrativa, **no un cálculo de campo tetragonal ni un Hamiltoniano atomístico validado**.

La restricción de red es `K_red*mean((r_j-R_red)**2)/2` por centro, y opcionalmente se penaliza deformar la esfera mediante `K_shape*sum((r_j-mean(r_j))**2)/2`. Ambas se reducen a la expresión radial original para seis radios iguales. Solo para dos Fe puede activarse `K_bridge*(mean(r_1j)+mean(r_2j)+L_bridge-R_bridge)**2/2`: es una **restricción de longitud de un fragmento idealizado**, no una descripción atomística de la red de Hofmann. Con `K_bridge=0` los centros son independientes. No hay intercambio magnético.

`protocol="nominal"` emplea `R_red(T)=R_0K+(R_300K-R_0K)*T/300` y barridos desde `nominal_start` hasta `nominal_end` en pasos `delta_T`. Aquí **T es únicamente una etiqueta geométrica**, nunca se usa `kBT`. `protocol="distance"` barre directamente desde `R_start` hasta `R_end` en pasos `delta_R`; la columna `nominal_T_K` será vacía. La secuencia de calentamiento recorre el camino inverso conservando el estado final del enfriamiento.

En **cada geometría fija** se evalúan las dos ramas electrónicas y se asigna `HS` si `E_HS <= E_LS`, o `LS` en caso contrario; `E_ground=min(E_HS,E_LS)` por centro. **El espín nunca se propone aleatoriamente.** En cada intento MC se copia la geometría, se desplaza temporalmente un único radio elegido al azar, se recalcula el estado fundamental y se acepta solo si la **energía total** propuesta es menor o igual; si se rechaza no cambian coordenadas ni espín. Al cambiar la red, se recalcula el estado de la misma geometría antes del primer intento. **Un barrido MC = `N_Fe` intentos (radial) o `6*N_Fe` intentos (seis radios)**; esto no equivale a segundos. La referencia cuasiestática minimiza por separado cada configuración HS/LS y compara las energías optimizadas; en el modo de seis ligandos se usan varios puntos de partida numéricos, sin garantía matemática de mínimo global. `sweeps_fast`, `sweeps_slow`, `replicas` y `seed` son configurables; para ensayar condiciones geométricas iniciales diferentes se puede ejecutar otra vez con `initial_branch="LS"` en vez de `"HS"` (solo cambia el **radio inicial**, no fuerza el espín). Los tres ejemplos utilizan `P_LS=3.9 eV` como valor **ilustrativo** que permite alcanzar el cruce: no es un dato experimental. Para investigar el atrapamiento, probar `P_LS=4.7 eV` en otro barrido y comparar contra el equilibrio; no cambiar el algoritmo para forzar la transición.

## Archivos generados

- `config_used.json`: configuración efectiva completa.
- `equilibrium.csv`: estado fundamental cuasiestático por distancia externa.
- `spin_surfaces.csv`: mínimos de HS/LS (o las cuatro combinaciones para dos Fe) por distancia.
- `trajectories.csv`: cada réplica y paso con distancias medias `r1_A`, `r2_A`, `delta1_eV`, `delta2_eV` y estados `spin1_HS`, `spin2_HS` (HS=1, LS=0); en modo de seis radios también `fe1_n1_A` … `fe1_n6_A` y sus coordenadas relativas `fe1_n1_x_A`, `fe1_n1_y_A`, `fe1_n1_z_A`, etc. (también en `equilibrium.csv` y `spin_surfaces.csv`).
- `attempts.csv`: **cada intento** con centro y ligando elegido, distancia antes/propuesta/después, energías antes/propuesta/después, estado electrónico antes/propuesto/después y decisión de aceptación. Puede ocupar decenas de MB en barridos largos; no contiene un tiempo físico.
- `averages.csv`: promedios entre réplicas; un espín medio fraccionario no significa un estado electrónico intermedio.
- `convergence.csv`: energía durante los barridos MC en cada paso.
- `surfaces.png` (incluye `Δ_o(r)`), `trajectories.png`, `convergence.png`, en modo nominal `nominal_distance.png` y en modo de seis radios `ligand_distances.png`.

**Interpretación:** el cruce a geometría fija, el cruce entre energías *ya minimizadas* y el cambio observado durante MC pueden ocurrir en diferentes puntos. Los pequeños movimientos de descenso pueden quedar bloqueados ante una barrera: **que la referencia cambie a LS no garantiza que el MC geométrico alcance LS**. Un salto de geometría en el equilibrio no requiere ruido térmico; un bucle producido por Monte Carlo a 0 K puede ser atrapamiento algorítmico, **no una predicción de histéresis experimental en K/min**. La persistencia o desaparición del bucle con más barridos debe comprobarse. La aproximación de estado electrónico instantáneo tampoco calcula barreras o tasas físicas de interconversión entre multiplicidades. No se simulan vibraciones térmicas, susceptibilidad ni momento magnético.

El cruce `r_c` de `surfaces.png` se calcula **solo para seis radios iguales**: con radios distintos, la condición de cruce depende de los seis valores y no hay una única distancia crítica escalar. Las gráficas no representan niveles multiplete completos de Tanabe–Sugano.
