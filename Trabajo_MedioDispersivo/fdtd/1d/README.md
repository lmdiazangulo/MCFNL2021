# Simulation of a dispersive material

**Alumnos:** Ignacio Baena Jiménez y Gregorio García Valladares

## Modelo

Usamos el método FDTD para resolver un  medio dispersivo caracterizado por una permitividad que depende de la frecuencia. La permitividad suele ser una función dada por una suma de pares de polo de Debye o de Lorentz. 

En nuestro caso, va a ser la suma de pares de polos complejos conjugados
$$
\varepsilon(\omega) = \varepsilon_0\varepsilon_\infty + \varepsilon_0 \sum_{p=1}^P c_p / (j\omega-a_p) +c_p^*/(j\omega-a_p^*).
$$
Esto nos permite definir unas ecuaciones auxiliares, que dependen de unas nuevas variables $\vec{J}_p(\omega)$ y $\vec{J}'_p(\omega)$, que en el dominio del tiempo son
$$
\frac{d}{dt}\vec{J}_p(t) -a_p \vec{J}_p(t) = \varepsilon_0 c_p \frac{d}{dt}\vec{E}(t) \\
\frac{d}{dt}\vec{J}'_p(t) -a_p^* \vec{J}'_p(t) = \varepsilon_0 c_p^* \frac{d}{dt}\vec{E}(t).
$$
Dado que $\vec{E}(t)$ es real, $\vec{J}'_p(t)=\vec{J}_p^*(t)$, con lo que solo necesitamos trabajar con una de las dos ecuaciones diferenciales anteriores. Esta ecuación, junto a las ecuaciones rotacionales habituales de Maxwell, nos permiten implementar un FDTD estándar usando el algoritmo de Yee: 
	1) Actualizo el campo eléctrico en $n\Delta t$, que depende del campo eléctrico y las corrientes en el instante $(n-1)\Delta t$ y del campo magnético en  $(n-1/2)\Delta t$.
	2) Actualizo las corrientes de polarización a $n\Delta t$, que depende de su valor anteriores y de los campos eléctricos.
	3) Actualizo el campo magnético a $(n+1/2)\Delta t$ con Yee estándar.
	
## Archivos

Tenemos un archivo principal **mainDisp.py** que va a llamar a otros archivos para crear el mallado **mesh.py**, resolver las ecuaciones y graficarlas **viewerDisp.py**.  La resolución se hace mediante distintas clases, dependiendo de la zona dispersiva que queramos:
1) No dispersivo: **solver.py**
2) Toda la región dispersiva: **DispMedia.py**
3) Región dispersiva: **PartialDisp.py**

Además, se añaden archivos .gif donde se recogen animaciones para distintos casos. 

