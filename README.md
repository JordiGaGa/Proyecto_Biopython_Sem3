# Proyecto_Biopython_Sem3 : Efectos en la regulación génica del estrés inducido a Escherichia coli por la restricción de Mg2+

Este es un proyecto colaborativo el cual contiene distintos scripts de Python para el análisis de los datos de expresión de Escherichia coli bajo condiciones de estrés (en presencia de baja cantidad de magnesio) para mas información consulte: [Proyecto_Biopython_Sem3/docs/Reporte_Proyecto.md](https://github.com/JordiGaGa/Proyecto_Biopython_Sem3/blob/main/docs/Reporte_Proyecto.md).

## Uso 

Para empezar a trabajar se deben seguir los siguientes pasos:

```bash
# Crear un ambiente de conda en donde se instalen todas las librerias necesarias
conda env create -f environment.yml

# Para activar el ambiente...
conda activate Proyecto_Bio  

# Para desactivarlo
conda deactivate
```

Se puede correr dentro de línea de comando el conjunto de scripts o llamar de manera individual cada una de las funciones que lo componen.

### Desde línea de comando

Para correr el programa desde línea de comando, se toma como argumentos: el path hacia el archivo que contiene los datos de expresión, el nombre de las columnas que contienen los datos de control, el nombre de las condiciones a probar, el correo del usuario y el organismo del cual se obtuvieron los datos, así como otros parámetros opcionales, tal y como se muestra a continuación:

```bash
python main.py -i 'archivo' -cn [col_control1,col_controln] -st [col_state1,col_staten] -co 'correo' -org 'Organismo' -pv double -d 'directorio/' -pl 'plot' -t 'x' -sv
```

Donde: 

Los parámetros obligatorios son:
- `'archivo'` es el nombre del archivo que contiene los valores de expresión. El archivo debe contener los valores de un gen por línea.
- `[col_control1,col_control2]` son los nombres de las columnas control, separados por coma, en formato de lista.
- `[col_state1,col_state2]` son los nombres de las columnas, separados por coma, con las condiciones que se van a comparar contra el control.
- `'correo'` es el correo del usuario que se va a utilizar para realizar la búsqueda en Entrez.
- `'organismo'` es el organismo al cual pertenecen los datos del archivo de expresión.

Los parámetros opcionales son:
- `double` es el valor (del p-value ajustado) que se va a utilizar como corte para la selección de genes sobre y sub expresados.
- `'directorio/'` es la ruta al directorio en donde se van a almacenar los archivos de salida.
- `'plot'` indica el tipo de grafica a realizar con los resultados ya sea una de pie ('pie'), una gráfica para mostrar la distribución ('dist'), ambas ('all') o ninguna (None).
- `'x'` si se hace un gráfico de distribución indica si este es un bloxplot ('box') o un violinplot ('vl').
-  El último parámetro `-sv` permite al usuario decidir si quiere que se guarden tanto los archivos como las imágenes. Si este parámetro está activo pero no se ingresa un directorio se guardaran todos los archivos en la ubicación desde donde se corra el script.

### En python

Para correr la función `analisis_diferencial` del archivo `Expresion_diferencial.py`, como una función de Python, se toman como argumentos: el path hacia el archivo que contiene los datos de expresión, además de un diccionario que el nombre de las columnas que contienen los datos de control y las condiciones a probar. Ejemplo de uso:

```python
analisis_diferencial(table_csv:'./path_al_archivo_.csv', samples:{'control':[columna_control_1,columna_control_2], 'states': [conlumna_condicion_1, columna_condicion_2]})
```

Para correr la función `id_Gene` del archivo `Informacion_Genes.py`, como una función de Python, se toman como argumentos: El nombre de un gen, el organismo al que pertenece  junto con el correo que se va a utilizar para la búsqueda en entrez . Ejemplo de uso:

```python
Id_Gene(gen = 'yecR',org = 'Escherichia coli',co ='ejemplo_correo@unam.mx')
```

Para correr la función `gen_fuction_tag` del archivo `Informacion_Genes.py`, como una función de Python, se toman como argumentos: Un id especifico perteneciente a un gen, junto con el correo que se va a utilizar para la búsqueda en entrez . Ejemplo de uso:

```python
gen_function_tag(id = 946386, co = 'ejemplo_correo@unam.mx')
```

Para correr la función `pie_expresion_plot` del archivo `Data_viz.py`, como una función de Python, se toman como argumentos: la serie que obtienes como resultado de usar la función `analisis_diferencial` del archivo `Expresion_diferencial.py`, si se desea guardar la imagen generada y el directorio en donde se almacenará dicha imagen.

```python
pie_expresion_plot(df= pd.Series, save = True, output_dir: '.')
```

Para correr la función `expresion_dist_plot` del archivo `Data_viz.py`, como una función de Python, se toman como argumentos: la serie que obtienes como resultado de usar la función `analisis_diferencial` del archivo `Expresion_diferencial.py`, el tipo de gráfico que se quiere generar 'box' si es un box_plot o 'vl' si es un violin_plot, además tambien tiene la opción si se desea guardar la imagen generada y el directorio en donde se almacenará dicha imagen.

```python
expresion_dist_plot(df: pd.Series, graph: 'box', save : True, output_dir: '.'):
```
## Salida

El script main tiene varias opciones de corrida que pueden dar diferentes salidas dependiendo las necesidades del usuario, a continuacion se explicaran algunas de ellas

### Línea de comando
#### Default
 ```bash
 python main.py -i 'archivo' -cn col_control1,col_controln -st col_state1,col_staten -co 'correo' -org 'Organismo'
 ```
La salida de esta línea de comando, es visualizar a pantallas las dos gráficas default que tiene el programa 'pie' y 'box'
#### Guardar todo
 ```bash
 python main.py -i 'archivo' -cn col_control1,col_controln -st col_state1,col_staten -co 'correo' -org 'Organismo' -sv -d 'dir/'
 ```
Incorporar el parámetro 'sv' a la línea de comandos, indica al programa que se van a salvar todos los archivos generados en el programa, tantos archivos csv con datos de expresión diferencial y jpg con las gráficas generadas. Por otro lado el parámetro 'd' es la ruta al directorio donde guardar los archivos siendo predeterminadamente el directorio donde se encuentra el usuario

#### Opciones gráficas
Por defecto se imprimen y/o generan dos gráficas, una gráfica de pie con los porcentajes de los genes sobre y sub expresados recuperados; además de una gráfica de distribución ya sea boxplot o violin plot. Sin embargo, el código le da la libertad al usuario a decidir que gráficas generar con el parámetro 'pl'.  


 ```bash
 python main.py -i 'archivo' -cn col_control1,col_controln -st col_state1,col_staten -co 'correo' -org 'Organismo' -pl 'pie'
```
La salida de esta línea es visualizar sólo la gráfica de pie


```bash
 python main.py -i 'archivo' -cn col_control1,col_controln -st col_state1,col_staten -co 'correo' -org 'Organismo' -pl 'dist' -t 'box'
``` 
```bash
 python main.py -i 'archivo' -cn col_control1,col_controln -st col_state1,col_staten -co 'correo' -org 'Organismo' -pl 'dist' -t 'vl'
```
Para el caso de la gráfica de distribución, el usuario debe decidir cuál gráfica generar con el parámetro 't' siendo 'box' por defecto generando un boxplot, pudiendo ser cambiado a violin indicándolo con 'vl'

#### Elección de p_value
```bash
python main.py -i 'archivo' -cn col_control1,col_controln -st col_state1,col_staten -co 'correo' -org 'Organismo' -pv 0.05
```
El usuario puede escoger qué valor de corte tomar para el p_value ajustado y filtrar los datos siendo por defecto 0.005

### Python

- La función `analisis_diferencial` del archivo `Expresion_diferencial.py` regresa un dataframe con los valores de expresión diferencial (sólo aquellos que muestren un p-value ajustado menor o igual 0.005 o al introducido por el usuario). 
- La función `id_Gene` del archivo `Informacion_Genes.py` regresa el id de la base de datos de gene del gen introducido.
- La función `gen_fuction_tag` del archivo `Informacion_Genes.py` regresa un diccionario que contiene como llaves: El id del gen, su nombre, la descripción de la proteína (en caso de tenerla) y el locus tag, además de entregar el valor correspondiente a dicha llave.
- La función `pie_expresion_plot` del archivo `Data_viz.py` regresa un gráfico de pie con el porcentaje de genes diferencialmente expresados que se sobreexpresan y subexpresan.
- La función `expresion_dist_plot` del archivo `Data_viz.py` regresa un box-plot o un violin-plot con la distribución de los genes diferencialmente expresados que se sobreexpresan y subexpresan.
## Datos 

El script está diseñado para operar con los datos del análisis de expresión extraídos de la base de dato GEO *[Series GSE276379](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276379)* que contiene el perfil de expresión de E. coli en condiciones de estrés bajo una pequeña cantidad de magnesio.

## Metadatos y documentación  

Este README ofrece información de uso básico de los scripts, tanto del modo de uso, así como del formato de salida. Para obtener información más detallada acerca de los datos utilizados puedes acceder al repositorio [Small-proteins-induced-under-low-magnesium-stress](https://github.com/yadavalli-lab/Small-proteins-induced-under-low-magnesium-stress).

## Código fuente

El código fuente está disponible en este repositorio [Proyecto_Biopython_Sem3/src at main · JordiGaGa/Proyecto_Biopython_Sem3](https://github.com/JordiGaGa/Proyecto_Biopython_Sem3/tree/main/src). Se acoge con satisfacción cualquier contribución o sugerencia a través de solicitudes pull request.

## Términos de uso

Este script está disponible bajo la licencia [MIT License](https://github.com/JordiGaGa/Proyecto_Biopython_Sem3/blob/main/LICENSE). Consulte el archivo LICENSE para obtener más detalles.

## Cómo citar

Si utiliza este script en su trabajo, por favor cite: Jordi et al. (2024). *Proyecto_Biopython_Sem3* (Versión 2.7.5) [Repositorio GitHub]. GitHub. https://github.com/JordiGaGa/Proyecto_Biopython_Sem3


## Contáctenos

Si tiene problemas o preguntas, por favor abra un issue en este repositorio o póngase en contacto con nosotros en:
- Armando Gael Gónzalez Trapaga (aggonzal@lcg.unam.mx)
- Jordi García Garcés (jordigg@lcg.unam.mx)
- Jocelyn Trujillo Gutierrez (jocelynt@lcg.unam.mx)
