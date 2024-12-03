# Conclusiones 
Durante el desarrollo de este proyecto aprendimos a llevar una organización para un proyecto colaborativo en git hub, teniendo en cuenta el uso de issues, tags y ramas externas al main.

Otro punto muy importante que se aprendió en este proyecto es el de cómo acceder a campos más específicos de la base de datos "gene" para poder extraer los elementos recabados en el dataframe como lo son el ID del gen, la descripción de la proteína o el locus tag.

# ¿Cómo se podría mejorar el proyecto? 
Algunos de los puntos que podrían hacer más óptimo o enriquecedor el proyecto son los siguientes:

1. Hacer que las funciones de el script "Informacion_Genes" puedan recibir directamente la lista con los genes a los cuales se les va a sacar su ID para el caso de la función "Id_Gene", o el diccionario con el nombre del gen y su ID para la función "gen_function_tag". Esto con el fin de ser más amigable con el usuario en lugar de que se le vayan pasando los argumentos desde el main a través de una comprehension list, lo cual haría que pudieran llamarse directamente a las funciones sin tener que hacer uno su propia comprehension list en caso de no querer correr todo el main.
2. Podría más adelante implementarse una nueva función la cual realice un análisis de enriquecimeinto funcional de los genes subexpresados y sobreexpresados que sean recabados en los procesos previos, esto con el fin de poder caracterizar a mayor profundidad qué vías o procesos podrían estar viéndose afectados por las condiciones analizadas en el experimento realizado.   

# ¿Cómo se gestionó el trabajo en equipo? 
Se divideron las tareas de manera tal que cada quien se encargaría de al menos un script de los realizados, asimismo, todos aportaron a la documentación tanto en el archivo Reporte_Proyecto.md como en el README.md explicando lo que sus scripts hacían, cómo correr los códigos o poniendo características de los archivos almacenados en la carpeta data.