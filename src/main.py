# ===========================================================================
# =                            imports
# ===========================================================================
import argparse
import os
import pandas as pd 
from Expresion_diferencial import analisis_diferencial
from Informacion_Genes import Id_Gene, gen_function_tag
from logger import getlogger
from Utils import outliers_dif_exp
from Data_viz import pie_expresion_plot, expresion_dist_plot
# ===========================================================================
# =                            functions
# ===========================================================================
def process_list(argline):
    if os.path.isfile(argline):
        with open(argline) as file:
            return [ int(txid) for line in file.readlines() if (txid:=line.strip())]
    else:
        return [str(name) for name in argline.split(',') if name]
# ===========================================================================
# =                            main
# ===========================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Recibe una matrix de conteos y regresa las funciones de los genes mas y menos expresados')
    
    parser.add_argument('-i','--input', type=str, required=True,
                        help='Ruta a la matriz de conteo') # '../data/GSE276379_RNASeq_kallisto.csv'
    parser.add_argument('-cn','--control',type=process_list, required=True,
                        help='Lista con el nombre de las columnas que son el control en el dataframe') # n1,n2
    parser.add_argument('-st','--states',type=process_list, required=True,
                        help='Lista con el nombre de las columnas que son el estado alternativo en el dataframe') # low-mg1,low-mg2
    parser.add_argument('-sv','--save', action='store_true',
                        help='Salvar o no el dataframe final')
    parser.add_argument('-co', '--correo', type=str, required=True,
                        help='Correo electronico para acceder a Entrez') # jordigg@lcg.unam.mx 
    parser.add_argument('-org', '--organismo', type=str, required=True,
                        help='Nombre del organismo del cual se hizo el analisis de expresion diferencial') # 'escherichia coli'
    parser.add_argument('-pv','--adjusted_p_value',type=float, default=0.005,
                        help='Valor del p value ajustado para filtrar los datos')
    parser.add_argument('-pl','--plot', type=str,choices=['all','pie','dist','None'],default='all',
                        help='Opciones de plots')
    parser.add_argument('-d','--dir_output',default='.',
                        help='Directorio donde guardar el archivo con todos los datos de expresion')

    args, _ = parser.parse_known_args()

    
   
    if args.plot == 'all' or args.plot == 'dist':
        parser.add_argument('-t','--type_dist',choices=['box','vl'],default='box',
                            help='Tipo de grafica de distribucion')
    
    args = parser.parse_args()
    
    logger = getlogger()
    
    try:
        logger.debug('Creando diccionario')
        samples = {'states':args.states, 'control':args.control}
        
        logger.info('Ejecutando script de analisis diferencial')
        df = analisis_diferencial(table=args.input,samples=samples, p_value=args.adjusted_p_value)
        #Añadir Resto de Procesos
        
    except Exception as e:
        raise(f'Error{e}')
    #Extraemos los outliaers 
    logger.info('Extrayendo genes outlaiers')
    genes_out = outliers_dif_exp(df)
   
    # Llamada a funcion Id_Gene que regresa los IDs de los genes a buscar
    
    logger.info('Empezando obtencion de IDs de genes')
    gen_id_mas = {gen:Id_Gene(gen, args.organismo, args.correo) for gen in genes_out['Sobreexpresado']}
    gen_id_menos = {gen:Id_Gene(gen, args.organismo, args.correo) for gen in genes_out['Subexpresado']}

    # Llamada a funcion gen_function_tag que crea los dataframes con ID, gen, descripcion de la proteína y locus tag
    
    logger.info('Empezando creacion de dataframes')
    #dataframe_mas = pd.DataFrame.from_dict([ig.gen_function_tag(i, args.correo) for i in gen_id_mas.values()])
    dataframe_mas = pd.DataFrame.from_dict([res for res in (gen_function_tag(i, args.correo) for i in gen_id_mas.values()) if res is not None])

    #dataframe_menos = pd.DataFrame.from_dict([ig.gen_function_tag(i, args.correo) for i in gen_id_menos.values()])
    dataframe_menos = pd.DataFrame.from_dict([res for res in (gen_function_tag(i, args.correo) for i in gen_id_menos.values()) if res is not None])

    logger.info('DataFrames listos')
    # Creacion de los archivos dataframe a csv
    if args.save and not df.empty:
        
        df.to_csv(os.path.join(args.dir_output,'full_data.csv'))        
        dataframe_mas.to_csv(os.path.join(args.dir_output,'Genes_Sobreexpresados.csv'), index=False)  
        dataframe_menos.to_csv(os.path.join(args.dir_output,'Genes_Subexpresados.csv'), index=False) 

    
    logger.info('Comenzando plot')
    
    if args.plot != 'None':
        if args.plot == 'all' or args.plot == 'pie':
            pie_expresion_plot(df,save=args.save,output_dir=args.dir_output)
        if args.plot == 'all' or args.plot == 'dist':
            expresion_dist_plot(df,args.type_dist,save=args.save,output_dir=args.dir_output)

    


