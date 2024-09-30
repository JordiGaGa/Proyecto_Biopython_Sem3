# [Efectos en la regulación génica del estrés inducido a Escherichia coli por la restricción de Mg2+]

Nombre:  Armando Gael Gónzalez Trapaga (aggonzal@lcg.unam.mx)   
Nombre:  Jordi García Garcés (<jordigg@lcg.unam.mx>)   
Nombre:   Jocelyn Trujillo Gutierrez  (<jocelynt@lcg.unam.mx>)  

Fecha:  10/09/2024
## Introducción

El magnesio es uno de los cationes mas abundantes dentro de  las células [referencia Romani AM]. Este ion tiene varios propósitos fundamentales como estabilizar complejos de macromoléculas, del mismo modo también puede actuar como cofactor en una variedad de complejos enzimáticos, así como influye en la regulación de la expresión genética. Enfocándonos en este último punto, el magnesio actúa en 5 puntos claves: unión al DNA y mantenimiento de su estabilidad (promoviendo la unión de factores de transcripción a las regiones promotoras, facilitando abrir su estructura en doble hélice para iniciar la transcripción), catálisis enzimática (cofactor en RNA-pol, varias ribozimas y adenyltransferasas las cuales ayudan a la maduración del tRNA), interacciones proteína-proteína (facilita la interacción entre las subunidades de la RNA-pol, también está presente en la formación de macro-complejos que regulen el DNA), regulación post-transcripcional (splicing y formación de los residuos de pseudouridina en rRNA y tRNA) y respuesta al estrés celular (modulando la actividad de ciertos factores de transcripción responsables de la respuesta a estrés, así como afectando la estabilidad de los mRNAs inducidos por estrés). [referencia Xue Li ].

Otros cationes pueden remplazar al Mg+2 en algunos de estos procesos , pero hay otros que están totalmente restringidos al uso de este catión [referencia Groisman], indicando 2 cosas: La primera es vías alternas para la realización de ciertos procesos, lo ayuda a la supervivencia en condiciones limitadas del ion y la otra es que la célula no puede vivir sin el ion, por ende limitar su presencia en el medio induciría el suficiente estrés para permitir cambios en la expresión genética y activar estas vías alternas.

## Planteamiento del problema

Se busca analizar la diferencia en los niveles de expresión de E.coli en condiciones normales y bajo estrés (expuestas a una concentración pequeña de magnesio), esto con el fin de identificar aquellos genes que se sobreexpresen y subexpresen bajo dicha condición, para posteriormente poder hacer una anotación funcional e identificar si existe alguna relación entre las funciones de los genes más y menos expresados.

## Planteamiento de Hipótesis

Aquellos genes que se sobreexpresen podrían estar ligados a rutas metabólicas alternas que no dependan del magnesio, mientras que aquellos que se muestren subexpresados podrían estar relacionados a genes que regulen o disminuyan la entrada de magnesio a la célula.

## Metodología

<!-- [Identificar y describir los diferentes datos de entrada con los que se cuenta, así como de dónde fueron descargados, el formato de los mismos, y las columnas con las que cuenta. Especificar si se utilizará un servidor en particular para trabajar, o herramientas para el desarrollo de la solución del análsis. Formular las preguntas biológicas que se busca resolver con el análisis de los datos para determinar las tareas a realizar por cada una de ellas.]


### A. Servidor y software

> Servidor: 

> Usuario: 

> Software: 


### B. Datos de Entrada 
-->

Los datos de entrada fueron descargados desde la base de datos de GEO en NCBI y se encuentran en la RUTA DE LA CARPETA.

```
|-- data
|   |-- GSE276379_RNASeq_kallisto.csv
|   |-- GSE276379_RNASeq_kallisto_sin_length.tsv
|   |-- SraRunTable.txt
```

<!-- 
#### Metadatos de la carpeta de datos
 
> Versión/Identificador del genoma:  NC_000913.3

> Fecha de descarga: dd/mm/aaaa

>| Archivo | Descripción  | Tipo |
|:--      |:--           |:--  |
| coli_genomic.fna  | Secuencia de nucleotidos de E. coli  | Formato FastA |
| coli.gff.   | Anotación del genoma de E. coli  | Formato gff |
| coli_protein.faa | Secuencia de aminoacidos de las proteinas de E. coli | formato FastA|
| flagella_genes.txt | Genes con función relacionada al flagello en E. coli | lista |
| directorio.txt. | Archivo con nombres de personas | lista |

-->

#### Formato de los archivos

<!-- 

- `coli_genomic.fna` : formato FastA


```
>NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome
AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTG
GTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGAC
AGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGT
```

Formato: 

> a. La primera línea es información de la secuencia. Primero viene el identificador del genoma.

> b. Después vienen varias líneas con la secuencia de nuclótidos del genoma completo.



- `coli.gff`: anotación de features en el genoma


El contenido del archivo es:

```
##gff-version 3
#!gff-spec-version 1.21
#!processor NCBI annotwriter
#!genome-build ASM584v2
#!genome-build-accession NCBI_Assembly:GCF_000005845.2
##sequence-region NC_000913.3 1 4641652
##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=511145

NC_000913.3     RefSeq  region  1       4641652 .       +       .       ID=NC_000913.3:1.>
NC_000913.3     RefSeq  gene    190     255     .       +       .       ID=gene-b0001;Dbx>
NC_000913.3     RefSeq  CDS     190     255     .       +       0       ID=cds-NP_414542.>
NC_000913.3     RefSeq  gene    337     2799    .       +       .       ID=gene-b0002;Dbx>
NC_000913.3     RefSeq  CDS     337     2799    .       +       0       ID=cds-NP_414543.>

```

Formato: 

> a. Es un formato gff tabular, es decir, cada dato es separado por tabulador.
> 
> b. Cada renglón en el formato gff es una elemento genético anotado en el genoma, que se le denomina `feature`, éstos features pueden ser genes, secuencias de inserción, promotores, sitios de regulación, todo aquello que este codificado en el DNA y ocupe una región en el genoma de  E. coli.

> c. Los atributos de cada columna par cada elemento genético son

>```
1. seqname. Nombre del cromosoma
2. source. Nombre del programa que generó ese elemento
3. feature. Tipo de elemento
4. start. Posición de inicio
5. end. Posición de final
6. score. Un valor de punto flotante
7. strand. La cadena (+ , - )
8. frame. Marco de lectura
9.  attribute. Pares tag-value, separados por coma, que proveen información adicional
```
--->

#### Preguntas de investigación
> ¿Existe un cambio en la expresión de los genes en condiciones bajas de magnesio?
 1. Obtención de los datos de expresión de E. coli en condiciones de bajo magnesio y grupo control.
 2. Modificar y reestrcturar los datos para poder manipular.
 3. Normalización y análisis de los datos de expresión con el programa PyDESeq2.
 4. Dado los datos obtenidos de PyDESeq2 identificar si hay cambios en la expresión. 

> ¿Qué genes están notablemente sobreexpresados y subexpresados?
 1. Identificar aquellos genes con cambios muy notables en sus niveles de expresión.
 2. Clasificar de acuerdo a si están sobreexpresados o subexpresados.

> ¿A qué funciones biológicas están asociadas dichos genes? 
 1. Realizar la anotación funcional de aquellos genes diferencialmente expresados con el módulo de Entrez. 
 2. Graficación de los datos. 



## Resultados
 

<!-- ### X. Pregunta 

Archivo(s):     

Algoritmo: 

1. 

Solución: Describir paso a paso la solución, incluyendo los comandos correspondientes

```bash

```

-->




## Análisis y Conclusiones

 <!-- Describir todo lo que descubriste en este análisis -->


## Referencias
<!-- Registrar todas las referencias consultadas. Se sugiere formato APA. Ejemplo:
 
 [1] Frederick R. Blattner et al., The Complete Genome Sequence of <i>Escherichia coli</i> K-12.Science277,1453-1462(1997).DOI:10.1126/science.277.5331.1453 -> EJEMPLO
 
  -->
  
 [1] Groisman EA, Hollands K, Kriner MA, Lee EJ, Park SY, Pontes MH. Bacterial Mg2+ homeostasis, transport, and virulence. Annu Rev Genet. 2013;47:625-46. doi: 10.1146/annurev-genet-051313-051025. Epub 2013 Sep 20. PMID: 24079267; PMCID: PMC4059682.
 
 [2] Romani AM, Scarpa A. Regulation of cellular magnesium. Front Biosci. 2000 Aug 1;5:D720-34. doi: 10.2741/romani. PMID: 10922296.
 
 [3] Xue Li, Xiaobai Zhang, Miaomiao Zhang, Xi Luo, Tingting Zhang, Xianjin Liu, Renfei Lu, Yiquan Zhang, Environmental magnesium ion affects global gene expression, motility, biofilm formation and virulence of Vibrio parahaemolyticus, Biofilm, Volume 7, 2024, 100194, ISSN 2590-2075, doi:10.1016/j.bioflm.2024.100194.
 
