# ⚡ Simulador Electromagnético 2D · FDTD (Yee)

Simulador interactivo en un solo archivo HTML que resuelve las ecuaciones de
Maxwell en 2D mediante el esquema de **Yee (FDTD)** en modo **TM**, con
fronteras **CPML** y varios presets didácticos (reflexión, guiado, cristal
fotónico).

Sin dependencias, sin build step, sin servidor. Abres el HTML y funciona.

---

## 🎯 ¿Qué hace?

Resuelve en el dominio del tiempo las ecuaciones de Maxwell para el modo
transversal-magnético (TM: `Ez`, `Hx`, `Hy`) sobre una malla 2D, y muestra en
tiempo real:

- El campo eléctrico `Ez(x, z, t)` como mapa de color divergente
  (azul ↔ negro ↔ naranja).
- Las trazas temporales de varias **sondas** repartidas por el dominio.
- HUD con tiempo físico normalizado y número de paso.

Con los presets incluidos puedes observar, sin escribir una línea:

| Preset | Fenómeno físico observable |
|---|---|
| **Vacío** | Onda cilíndrica 2D. Amplitud `∝ 1/√r` lejos de la fuente. |
| **Lámina dieléctrica** | Reflexión `R = ((1−n)/(1+n))² ≈ 0.111` en cada cara; rebotes internos. |
| **Guía de onda** | Paredes PEC arriba/abajo; propagación en `x`. Superposición de modos. |
| **Cristal fotónico** | Red cuadrada de varillas `ε_r = 11.4`. Dispersión múltiple y, a ciertas frecuencias, *band gap* (reflexión casi total). |

La fuente se puede configurar como **Ez blanda** (dipolo puntual), **Ez dura**
(fuente ideal) u **onda plana** (línea horizontal).

---

## 🧠 ¿Cómo funciona?

### Física: ecuaciones de Maxwell en modo TM

En 2D, con propagación en el plano `(x, z)` y campos independientes de `y`, se
separan dos polarizaciones. Aquí se resuelve la **TM** (o *E-parallel*),
cuyas únicas componentes no nulas son `Ez`, `Hx`, `Hy`:
∂Ez/∂t = (1/ε) [ ∂Hy/∂x − ∂Hx/∂z − σ Ez ]
∂Hx/∂t = −(1/μ) [ ∂Ez/∂z ]
∂Hy/∂t = (1/μ) [ ∂Ez/∂x ]

text

Se trabaja en **unidades normalizadas** (`c = ε₀ = μ₀ = 1`) para evitar
factores de `10⁸` y hacer las constantes tratables a mano. La impedancia del
vacío vale `η = 1`, y la velocidad de fase en un medio es `1/√(ε_r μ_r)`.

### Numérico: esquema de Yee

La clave del método es que las componentes de `E` y `H` se sitúan en
**posiciones escalonadas** de la malla (grid de Yee):
Ez(i, j) en (i·dx, j·dz)
Hx(i, j) en (i·dx, (j+½)·dz)
Hy(i, j) en ((i+½)·dx, j·dz)

text

Con esta colocación, las diferencias finitas centradas son de segundo orden
**sin interpolación**, y la actualización es de tipo *leapfrog*:
H^{n+½} = H^{n−½} − (dt/μ) ∇×E^n
E^{n+1} = E^n + (dt/ε) ∇×H^{n+½}

text

El campo magnético se actualiza medio paso por delante del eléctrico, lo que
da al esquema estabilidad y conservación de energía sin disipación numérica
espuria.

**Estabilidad (CFL):** en 2D el paso temporal está acotado por
dt ≤ 1 / (c · √(1/dx² + 1/dz²)) con dx = dz = 1 → dt ≤ 1/√2

text

El simulador usa `dt = S · dx / (c·√2)` con `S` (Courant) ajustable entre 0.1
y 0.7. El valor por defecto es `S = 0.5`, holgadamente dentro del régimen
estable.

### Fronteras: CPML

Un dominio finito sin tratamiento absorbe mal y produce reflexiones parásitas
que arruinan cualquier observación de ondas salientes. Se usa **CPML**
(Convolutional PML, Roden & Gedney 2000) con perfil de conductividad
polinómico:
σ(ρ) = σ_max · (ρ/L)^m m = 3, L = 20 celdas
σ_max = −(m+1) · ln(R₀) / (2 · η · L) R₀ ~ 1e−6

text

Implementación práctica: la derivada espacial en el interior del PML se
sustituye por `∂/∂x → ∂/∂x + ψ`, donde `ψ` es una variable auxiliar que
obedece
ψ^{n} = b(x) · ψ^{n−1} + c(x) · (∂/∂x)^n
b(x) = exp(−σ(x)·dt/ε₀), c(x) = b(x) − 1

text

Con `κ = 1` y `α = 0` (los parámetros extra de CPML general se omiten porque
no aportan nada al caso isótropo). El resultado son reflexiones del orden de
`−60 dB` a `−80 dB`, imperceptibles en el mapa de color.

### Materiales

Cada celda lleva un objeto `Medium` con `ε_r`, `μ_r`, `σ` y un flag `pec`.
Los coeficientes de actualización se precalculan una vez por modelo:
ca = (1 − σ·dt/(2ε)) / (1 + σ·dt/(2ε)) (evolución de Ez con pérdidas)
cb = (dt/ε) / (1 + σ·dt/(2ε))

text

En conductores perfectos (`pec: true`) el campo se fuerza a `Ez = 0` tras
cada paso.

`μ` se promedia horizontal o verticalmente en las posiciones escalonadas de
`Hx` y `Hy` respectivamente, para no introducir asimetrías numéricas.

### Fuente

Pulso gaussiano modulado:
s(t) = A · sin(2π f (t − t₀)) · exp( −((t − t₀)/τ)² )
t₀ = 3/f, τ = 1.2/f

text

Ocupa unas 3–4 oscilaciones, suficiente banda para excitar varios modos pero
suficientemente estrecho como para no llenar el dominio de energía antes de
que llegue el frente de onda. Se inyecta **después** de actualizar `E`, para
que no sea sobreescrito por el solver.

---

## 🏗️ Arquitectura del código

El proyecto sigue una **arquitectura hexagonal** (puertos y adaptadores),
heredera directa del simulador sísmico hermano. La idea es que el dominio no
sepa que existe un navegador, y la UI no sepa que existe un solver.
┌─────────────────────────────────────────────────────────────┐
│ DOMINIO (puro) │
│ Medium · Source · Probe · SimConfig · MODEL_PRESETS │
│ Sin dependencias, testeable en Node. │
└─────────────────────────────────────────────────────────────┘
▲
│ (usa)
┌─────────────────────────────────────────────────────────────┐
│ INFRAESTRUCTURA (adaptadores) │
│ FDTDEngine → solver numérico │
│ FieldVisualizer → pinta Ez(x,z) en canvas │
│ ProbeRenderer → pinta trazas de sondas │
│ EventBus → mediador pub/sub │
└─────────────────────────────────────────────────────────────┘
▲
│ (orquesta vía eventos)
┌─────────────────────────────────────────────────────────────┐
│ APLICACIÓN (casos de uso) │
│ StartSimUC · StepSimUC · ResetSimUC │
│ UpdateSourceUC · AddProbeUC · UpdateModelUC │
└─────────────────────────────────────────────────────────────┘
▲
│
┌─────────────────────────────────────────────────────────────┐
│ UI │
│ EMUI → suscrita al bus, traduce clics a casos de uso │
└─────────────────────────────────────────────────────────────┘

text

**Regla de dependencias:** las flechas apuntan hacia dentro. La UI emite
casos de uso, los casos de uso llaman al motor, el motor emite eventos por el
bus, y la UI se limita a pintar lo que llega. En ningún punto la UI toca
`engine.ez` directamente.

### Por qué así

- **Testeabilidad**: el dominio y el motor se pueden importar en Node sin
  canvas ni `window`.
- **Sustituibilidad**: cambiar el esquema numérico (FDTD → FDFD, o Yee →
  pseudo-espectral) no toca la UI ni los presets.
- **Extensibilidad**: añadir TFSF o dispersión es tocar `step()` y los
  materiales, nada más.

---

## 🚀 Cómo usarlo

1. Descarga `electromagnetismo.html`.
2. Ábrelo en cualquier navegador moderno (Chrome, Firefox, Safari, Edge).
3. Selecciona un modelo en el desplegable.
4. Pulsa **Iniciar**.

### Controles

| Control | Qué hace |
|---|---|
| **Modelo** | Cambia de preset (vacío, dieléctrico, guía, cristal). |
| **Fuente** | Puntual (Ez blanda/dura) u onda plana. |
| **Frecuencia** | `f·Δx/c` normalizada. Sube → más oscilaciones por celda. |
| **Amplitud** | Escala la inyección. |
| **Courant** | `S = c·dt·√2 / Δx`. Baja → más preciso, más lento. |
| **Pasos por frame** | Compromiso entre velocidad visual y fluidez. |
| **Añadir sonda** | Coloca una sonda en posición aleatoria del dominio. |

### Cosas que merece la pena probar

1. **Cristal fotónico**: sube `frecuencia` de 0.06 a 0.11. Verás cómo la
   transmisión cae en picado al entrar en el *band gap*.
2. **Lámina dieléctrica**: cambia la fuente a `planeWave`. Incide desde
   arriba y observa el pulso reflejado y transmitido separarse.
3. **Courant = 0.3** en el dieléctrico: pulsos más limpios, dispersión
   numérica apreciablemente menor.
4. **Añade 5–6 sondas** en línea horizontal para ver el retardo entre
   llegadas (velocidad de fase).

---

## 📁 Estructura del proyecto
.
├── electromagnetismo.html ← todo el simulador en un archivo
├── README.md ← este archivo
└── (opcional) docs/
├── derivacion.md ← deducción paso a paso del esquema
└── cpml.md ← derivación del perfil de σ y constantes

text

El simulador es un único archivo autocontenido. No hay bundler, no hay
`node_modules`, no hay transpilación. La sección `<script>` está dividida en
bloques comentados (`// ==== DOMINIO ====`, `// ==== INFRAESTRUCTURA ====`,
etc.) para navegarlo rápido.

---

## ⚠️ Limitaciones y honestidad

Esto es una herramienta **didáctica**. Se han omitido conscientemente cosas
que un solver de producción sí tiene:

- **Solo 2D TM.** No hay modo TE (`Hz, Ex, Ey`) implementado, aunque es una
  extensión trivial (~20 líneas).
- **Sin TFSF.** La onda plana del preset es una línea inyectada, no un frente
  unidireccional perfecto. Verás reflejos de la fuente hacia el PML superior.
- **Sin dispersión material.** No hay Drude, Lorentz ni Debye. Para
  metamateriales con `ε(ω) < 0` hay que añadir una ecuación diferencial
  auxiliar (ADE) por polo.
- **Sin anisotropía.** Los materiales son isótropos.
- **Fuente puntual ideal.** No hay antena con impedancia, ni línea de
  transmisión, ni puerto de excitación realista.
- **CPML con `κ = 1, α = 0`.** Suficiente para isótropo sin ángulos rasantes
  extremos, insuficiente para medios muy conductores.
- **Sin absorbente de banda ancha.** El pulso gaussiano modulado da un ancho
  de banda limitado. Para análisis espectral hay que barrer frecuencias.

---

## 🔬 Extensiones posibles

Por orden de dificultad creciente:

1. **Modo TE.** Duplicar el solver con `Hz, Ex, Ey`. 20 líneas.
2. **Vector de Poynting.** `Sx = −Ez·Hy`, `Sy = Ez·Hx`. Muy visual.
3. **TFSF.** Onda plana unidireccional real. Unas 60 líneas, pero hay que
   separar el dominio en zona total y zona dispersada.
4. **Dispersión Drude/Lorentz.** Dos variables extra por celda dispersiva.
5. **Medios anisotrópos.** `ε` pasa a tensor 3×3 por celda.
6. **Guardar/cargar modelos.** Serializar `materials` a JSON.
7. **Exportar trazas a CSV.** Botón de descarga en las sondas.

---

## 📚 Referencias

- **Yee, K. S.** (1966). *Numerical solution of initial boundary value
  problems involving Maxwell's equations in isotropic media.* IEEE Trans.
  Antennas Propag. 14(3), 302–307.
- **Taflove, A. & Hagness, S. C.** (2005). *Computational Electrodynamics:
  The Finite-Difference Time-Domain Method.* 3rd ed. Artech House.
- **Roden, J. A. & Gedney, S. D.** (2000). *Convolutional PML (CPML): An
  efficient FDTD implementation of the CFS-PML for arbitrary media.*
  Microwave and Optical Technology Letters, 27(5), 334–339.
- **Bérenger, J.-P.** (1994). *A perfectly matched layer for the absorption
  of electromagnetic waves.* J. Comput. Phys. 114(2), 185–200.

---

## 📜 Licencia

agp3. El código está pensado para enseñar; úsalo,
modifícalo, rómpelo.

---

## 🔗 Proyectos relacionados

- **`simulador_sismico3.html`** — simulador elastodinámico 2D con el mismo
  esquema (Virieux, grid escalonado, leapfrog). Estructura del código y
  estilo de UI idénticos; sirve como comparación directa entre FDTD
  electromagnético y FD elastodinámico.

Ambos comparten la misma arquitectura hexagonal y el mismo bus de eventos, lo
cual hace que moverse de uno a otro sea casi inmediato
