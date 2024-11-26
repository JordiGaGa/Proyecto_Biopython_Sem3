# Proyecto_Biopython_Sem3 : Efectos en la regulación genica del estrés inducido a Escherichia coli por la restricción de Mg2+

Este es un proyecto colaborativo el cual contiene distintos scripts de Python para el análisis de los datos de expresión de Escherichia Coli bajo condiciones de estrés (en presencia de baja cantidad de magnesio) para mas información consulte: [Proyecto_Biopython_Sem3/docs/Reporte_Proyecto.md](https://github.com/JordiGaGa/Proyecto_Biopython_Sem3/blob/main/docs/Reporte_Proyecto.md).

## Uso (TODO)

Para empezar a trabajar se deben seguir los siguientes pasos:

```bash
# RECOMENDACION: Crear un ambiente de conda en donde instales todas las librerias necesarias
conda create --name p_biopython
conda activate p_biopython 
# OBLIGATORIO: Descargar las librerias necesarias para el uso de los proyectos
pip install pandas
pip install pydeseq2
pip install Biopython
```

Cada uno de los scripts acepta distintos argumentos.

El script `Expresion_diferencial.py`, toma como argumento el path hacia el archivo que contiene los datos de expresión, el nombre de las columnas que contienen los datos de control y el nombre de las condiciones a probar:

```
python Expresion_diferencial.py -i [archivo] -cn [col_control1,col_controln] -st [col_state1,col_staten] 
```

donde: 
- `[archivo]` es el nombre del archivo que contiene los valores de expresión. El archivo debe contener los valores de un gen por línea.
- `[col_control1,col_controln]` es el nombre de las columnas control en formato de lista.
- `[col_state1,col_staten]` es el nombre de las columnas con las condiciones que se van a comparar contra el control.

## Salida (TODO)

El script imprimirá la suma de los números en la consola.

## Pruebas (TODO)

El script incluye un conjunto de pruebas unitarias. Puede ejecutar estas pruebas con:

```
python -m unittest test_sumador.py
```

## Datos 

El script está diseñado para operar con los datos del análisis de expresión extraídos de la base de dato GEO *[Series GSE276379](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276379)* que contiene el perfil de expresion de E. coli en condiciones de estrés bajo una pequeña cantidad de magnesio.

## Metadatos y documentación  

Este README ofrece información de uso básico de los scripts, tanto del modo de uso, así como del formato de salida. Para obtener información mas detallada acerca sobre los datos utilizados puedes acceder al repositorio [Small-proteins-induced-under-low-magnesium-stress](https://github.com/yadavalli-lab/Small-proteins-induced-under-low-magnesium-stress).

## Código fuente (TODO)

El código fuente está disponible en este repositorio [Proyecto_Biopython_Sem3/src at main · JordiGaGa/Proyecto_Biopython_Sem3](https://github.com/JordiGaGa/Proyecto_Biopython_Sem3/tree/main/src). Se acoge con satisfacción cualquier contribución o sugerencia a través de solicitudes pull request.

## Términos de uso

Este script está disponible bajo la licencia [MIT License](https://github.com/JordiGaGa/Proyecto_Biopython_Sem3/blob/main/LICENSE). Consulte el archivo LICENSE para obtener más detalles.

## Como citar

Si utiliza este script en su trabajo, por favor cite: Jordi et al. (2024). *Proyecto_Biopython_Sem3* (Versión 1) [Repositorio GitHub]. GitHub. https://github.com/JordiGaGa/Proyecto_Biopython_Sem3


## Contáctenos

Si tiene problemas o preguntas, por favor abra un problema en este repositorio o póngase en contacto con nosotros en:
- Armando Gael Gónzalez Trapaga (aggonzal@lcg.unam.mx)
- Jordi García Garcés (jordigg@lcg.unam.mx)
- Jocelyn Trujillo Gutierrez (jocelynt@lcg.unam.mx)
