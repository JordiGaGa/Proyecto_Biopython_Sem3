# [Efectos en la regulación génica del estrés inducido a Escherichia coli por la restricción de Mg2+]

Nombre:  Armando Gael Gónzalez Trapaga (aggonzal@lcg.unam.mx)   
Nombre:  Jordi García Garcés (<jordigg@lcg.unam.mx>)   
Nombre:   Jocelyn Trujillo Gutierrez  (<jocelynt@lcg.unam.mx>)  

Fecha:  10/09/2024

## Introducción

El magnesio es uno de los cationes mas abundantes dentro de  las células [2]. Este ion tiene varios propósitos fundamentales como estabilizar complejos de macromoléculas, del mismo modo también puede actuar como cofactor en una variedad de complejos enzimáticos, así como influye en la regulación de la expresión genética. Enfocándonos en este último punto, el magnesio actúa en 5 puntos claves: unión al DNA y mantenimiento de su estabilidad (promoviendo la unión de factores de transcripción a las regiones promotoras, facilitando abrir su estructura en doble hélice para iniciar la transcripción), catálisis enzimática (cofactor en RNA-pol, varias ribozimas y adenyltransferasas las cuales ayudan a la maduración del tRNA), interacciones proteína-proteína (facilita la interacción entre las subunidades de la RNA-pol, también está presente en la formación de macro-complejos que regulen el DNA), regulación post-transcripcional (splicing y formación de los residuos de pseudouridina en rRNA y tRNA) y respuesta al estrés celular (modulando la actividad de ciertos factores de transcripción responsables de la respuesta a estrés, así como afectando la estabilidad de los mRNAs inducidos por estrés). [3].

Otros cationes pueden remplazar al Mg+2 en algunos de estos procesos , pero hay otros que están totalmente restringidos al uso de este catión [1], indicando 2 cosas: La primera es vías alternas para la realización de ciertos procesos, lo ayuda a la supervivencia en condiciones limitadas del ion y la otra es que la célula no puede vivir sin el ion, por ende limitar su presencia en el medio induciría el suficiente estrés para permitir cambios en la expresión genética y activar estas vías alternas.

## Planteamiento del problema

Se busca analizar la diferencia en los niveles de expresión de E.coli en condiciones normales y bajo estrés (expuestas a una concentración pequeña de magnesio), esto con el fin de identificar aquellos genes que se sobreexpresen y subexpresen bajo dicha condición, para posteriormente poder hacer una anotación funcional e identificar si existe alguna relación entre las funciones de los genes más y menos expresados.

## Planteamiento de Hipótesis

Aquellos genes que se sobreexpresen podrían estar ligados a rutas metabólicas alternas que no dependan del magnesio, mientras que aquellos que se muestren subexpresados podrían estar relacionados a genes que regulen o disminuyan la entrada de magnesio a la célula.

## Metodología

<!-- [Identificar y describir los diferentes datos de entrada con los que se cuenta, así como de dónde fueron descargados, el formato de los mismos, y las columnas con las que cuenta. Especificar si se utilizará un servidor en particular para trabajar, o herramientas para el desarrollo de la solución del análsis. Formular las preguntas biológicas que se busca resolver con el análisis de los datos para determinar las tareas a realizar por cada una de ellas.]
-->

### A. Servidor y software

> Servidor: 
- NA

> Usuario: 
-  NA

> Software: 
- Python 


### B. Datos de Entrada 


Los datos de entrada fueron descargados desde la base de datos de GEO en NCBI y se encuentran en la RUTA DE LA CARPETA.

```
|-- data
|   |-- GSE276379_RNASeq_kallisto.csv
|   |-- GSE276379_RNASeq_kallisto_sin_length.tsv
|   |-- SraRunTable.txt
```


#### Metadatos de la carpeta de datos
 
> Versión/GEO accession: **[GSE276379](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276379)**

> Organismo :  **_Escherichia_ _coli_**

> Fecha de descarga: 09/2024

| Archivo                                  | Descripción                                                                           | Tipo        |
| :--------------------------------------- | :------------------------------------------------------------------------------------ | :---------- |
| GSE276379_RNASeq_kallisto.csv            | Archivo con los datos de expresión en distintas condiciones y la longitud de cada gen | Formato csv |
| GSE276379_RNASeq_kallisto_sin_length.tsv | Archivo con los datos de expresión en distintas condiciones                           | Formato tsv |
| SraRunTable.txt                          | Tabla de metadatos sobre cada una de las condiciones                                  | Formato csv |

#### Formato de los archivos


- `GSE276379_RNASeq_kallisto_sin_length.tsv` : formato tsv

```
	low-mg1	low-mg2	ns1	ns2
thrL	39.8686	40.0939	56.2871	61.5181
thrA	26.5438	21.487	341.707	298.131
thrB	25.3756	19.9654	255.076	237.161
thrC	25.6537	20.6849	310.112	285.386
yaaX	17.8188	14.4786	32.6503	34.8904
yaaA	29.9852	24.5718	60.3333	61.0447
```

Formato: 

> a. La primera línea es información de las condiciones: 
- low-mg1: E. coli en medio con bajo magnesio. Replica 1
- low-mg2: E. coli en medio con bajo magnesio. Replica 2
- ns1: Grupo control. Replica 1
- ns2: Grupo control. Replica 2

> b. Después vienen varias líneas con cada gen y su nivel de expresión en cada condición.

- `SraRunTable.txt`: Metadatos de los análisis


El contenido del archivo es:

```
Run,Assay Type,AvgSpotLen,Bases,BioProject,BioSample,Bytes,cell_type,Center Name,Collection_Date,Consent,DATASTORE filetype,DATASTORE provider,DATASTORE region,Experiment,genotype,geo_loc_name_country,geo_loc_name_country_continent,geo_loc_name,Instrument,Library Name,LibraryLayout,LibrarySelection,LibrarySource,Organism,Platform,ReleaseDate,create_date,version,Sample Name,source_name,SRA Study,strain,treatment

SRR30561641,RNA-Seq,290,14043346823,PRJNA1156803,SAMN43501102,5466390979,Bacterial cell,"YADAVALLI, GENETICS, RUTGERS UNIVERSITY",missing,public,"fastq,run.zq,sra","gs,ncbi,s3","gs.us-east1,ncbi.public,s3.us-east-1",SRX25984781,Wild type,uncalculated,uncalculated,missing,Illumina HiSeq 4000,GSM8497425,PAIRED,cDNA,TRANSCRIPTOMIC,Escherichia coli,ILLUMINA,2024-09-05T00:00:00Z,2024-09-05T18:34:00Z,1,GSM8497425,Bacterial cell,SRP530775,K-12 MG1655,low magnesium
```

Formato: 

> a. La primera línea son los datos relacionados a cada una de las pruebas realizadas. (Aquellos que aparezcan como missing, uncalculated o public en los metadatos no serán descritos)
- Run: Identificador asociado a la corrida de secuenciación en el repositorio de datos de NCBI.
- Assay Type: Tipo de ensayo con el que obtuvieron los datos.
- AvgSpotLen: Longitud promedio de las lecturas 
- Bases: Número de Bases
- BIoProject: Identificador del proyecto biológico relacionado
- BioSample: Identificador de la muestra biológica 
- Bytes: Cantidad de datos cuantificados en bytes.
- cell type: Tipo de célula analizada
- Center Name: Nombre del laboratorio o centro que realizó el estudio
- Collection_Date: 
- Consent:
- DATASTORE filetype: Formato o tipo de archivo 
- DATASTORE provider: Plataforma para almacenar los datos
- DATASTORE region: Ubicación del centro de datos de los datos experimentales.
- Experiment: Identificador del experimento.
- genotype: Características del genotipo (wild type o modificado)
- geo_loc_name_country:
- geo_loc_name_country_continent:
- geo_loc_name: 
- Instrument: Modelo y nombre del secuenciador.
- Library Name: Nombra de la biblioteca de secuencación.
- LibraryLayout: Lecturas pareadas o solas.
- LibrarySelection: Selección de la biblioteca de datos en base a con lo que fue generada.
- LibrarySource: Biblioteca de donde se extrajeron los datos.
- Organism: Organismo o especie usada.
- Platform: Plataforma usada para el análisis
- ReleaseDate: Hora y fecha de liberación de los datos.
- create_date: Hora y fecha de registro en la base de datos.
- version: Número de versión
- Sample Name: Nombre de la muestra usada 
- source_name: Tipo de célula 
- SRA Study: EStudio registrado en Sequence Read Archive 
- strain: Cepa utilizada 
- treatment: Condición experimental bajo la que estuvo la muestra. 

> b. Después cada una de las filas está relacionada a las diferentes pruebas que se realizaron con sus características e información.

#### Preguntas de investigación

> ¿Existe un cambio en la expresión de los genes en condiciones bajas de magnesio?
1. Obtención de los datos de expresión de E. coli en condiciones de bajo magnesio y grupo control.
 2. Modificar y reestrcturar los datos para poder manipular.
 3. Normalización y análisis de los datos de expresión con el programa PyDESeq2.
 4. Dado los datos obtenidos de PyDESeq2 identificar si hay cambios en la expresión. 

> ¿Qué genes están notablemente sobreexpresados y subexpresados?
1. Identificar aquellos genes con cambios muy notables en sus niveles de expresión.
 2. Clasificar de acuerdo a si están sobreexpresados o subexpresados.
 3. Graficar.

> ¿A qué funciones biológicas están asociadas dichos genes? 
1. Buscar la descripción funcional de aquellos genes diferencialmente expresados con el módulo de Entrez. 

####Flujo del proyecto (Se tiene que actualizar)
![DiagramaFlujo](https://github.com/JordiGaGa/Proyecto_Biopython_Sem3/blob/main/docs/Diagrama_de_flujo.png)


## Resultados

### 1. ¿Existe un cambio en la expresión de los genes en condiciones bajas de magnesio? 

Archivo(s):     
	>`GSE276379_RNASeq_kallisto.csv`

Algoritmo: 

1.  Utilizar la libreria 

Solución: Describir paso a paso la solución, incluyendo los comandos correspondientes


```bash

```

Fig.1 pie_plot
### 2. ¿Qué genes están notablemente sobreexpresados y subexpresados? 

Archivo(s) (si se están corriendo como funciones individuales):     
	> full_data.csv

Algoritmo: 

1. Dado los resultados del analisis de expresion diferencial se separaran aquellos genes cuyo p-value haya sido significativo (0.005) en sobreexpresados y subexpresados.
2. De cada distribucion de valores se obtendran los outliers.
3. Se regresaran los outliers correspondientes para el siguente paso.

Solución: Describir paso a paso la solución, incluyendo los comandos correspondientes

```bash
import pandas as pd
from Informacion_Genes import Id_Gene
from Utils import outliers_dif_exp
from Data_viz import expresion_dist_plot
# Graficar la distribucion de la expresion de los genes diferencialmente expresados
expresion_dist_plot(pd.read_csv(full_data.csv)),'box')
# Separar y obtener los outliers de los genes mas sobre y sub expresados
genes_out = outliers_dif_exp(pd.read_csv(full_data.csv))
```

(Agregar figura)
Fig. 2 Box-plot de la distribución de los niveles de expresión de los genes 
### 3. ¿A qué funciones biológicas están asociadas dichos genes?

Archivo(s) (si se están corriendo como funciones individuales):     
	>

Algoritmo: 

1.  Dado los resultados de ...

Solución: Describir paso a paso la solución, incluyendo los comandos correspondientes


```bash

```

## Análisis y Conclusiones

De acuerdo con los resultados, **E. coli** adapta su mecanismo de regulación en respuesta a condiciones de baja disponibilidad de magnesio. Esta regulación diferencial muestra dos tendencias claras:
- **Genes subexpresados:** Predominantemente asociados con el **flagelo**. Esto podría indicar una reducción en la motilidad, probablemente como estrategia para conservar recursos energéticos en condiciones de estrés.
- **Genes sobreexpresados:** No presentan un patrón uniforme, lo que sugiere una activación de distintos mecanismos para adaptarse a la restricción de magnesio, como la regulación de transportadores, respuestas al estrés o ajuste del metabolismo.

## Referencias

 [1] Groisman EA, Hollands K, Kriner MA, Lee EJ, Park SY, Pontes MH. Bacterial Mg2+ homeostasis, transport, and virulence. Annu Rev Genet. 2013;47:625-46. doi: 10.1146/annurev-genet-051313-051025. Epub 2013 Sep 20. PMID: 24079267; PMCID: PMC4059682.
 
 [2] Romani AM, Scarpa A. Regulation of cellular magnesium. Front Biosci. 2000 Aug 1;5:D720-34. doi: 10.2741/romani. PMID: 10922296.
 
 [3] Xue Li, Xiaobai Zhang, Miaomiao Zhang, Xi Luo, Tingting Zhang, Xianjin Liu, Renfei Lu, Yiquan Zhang, Environmental magnesium ion affects global gene expression, motility, biofilm formation and virulence of Vibrio parahaemolyticus, Biofilm, Volume 7, 2024, 100194, ISSN 2590-2075, doi:10.1016/j.bioflm.2024.100194.
 
