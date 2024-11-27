# ===========================================================================
# =                            imports
# ===========================================================================
import argparse
import os
import pandas as pd 
from Expresion_diferencial import analisis_diferencial
import Informacion_Genes as ig
from logger import getlogger
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
    parser.add_argument('-sv_all','--save_full_data', action='store_true',
                        help='Salvar o no el dataframe final')
    parser.add_argument('-co', '--correo', type=str, required=True,
                        help='Correo electronico para acceder a Entrez') # jordigg@lcg.unam.mx 
    parser.add_argument('-org', '--organismo', type=str, required=True,
                            help='Nombre del organismo del cual se hizo el analisis de expresion diferencial') # 'escherichia coli'

    parser.add_argument('-p','--adjusted_p_value',type=float, default=0.005,
                        help='Valor del p value ajustado para filtrar los datos')

    args, _ = parser.parse_known_args()

    if args.save_full_data:
        parser.add_argument('-o_all','--output_full_data',default='./full_data.csv',
                            help='Directorio donde guardar el archivo con todos los datos de expresion')
        args = parser.parse_args()
    df = 0
    try:
        samples = {'states':args.states, 'control':args.control}
        
        df = analisis_diferencial(table=args.input,samples=samples p_value=arg.adjusted_p_value)
        #Añadir Resto de Procesos
       
        if args.save_full_data and not df.empty:
            df.to_csv(args.output_full_data)
        
    except Exception as e:
        raise(f'Error{e}')
    

    # Llamada a funcion Lista_genes que extrae los mas sobreexpresados y menos sobreexpresados 
    logger = getlogger('Informacion de genes')
    logger.info('Empezando creacion de lista de genes')
    genes_mas,genes_menos = ig.Lista_genes(df)
   
    # Llamada a funcion Id_Gene que regresa los IDs de los genes a buscar
    
    logger.info('Empezando obtencion de IDs de genes')
    gen_id_mas = {gen:ig.Id_Gene(gen, args.organismo, args.correo) for gen in genes_mas}
    gen_id_menos = {gen:ig.Id_Gene(gen, args.organismo, args.correo) for gen in genes_menos}

    # Llamada a funcion gen_function_tag que crea los dataframes con ID, gen, descripcion de la proteína y locus tag
    
    logger.info('Empezando creacion de dataframes')
    #dataframe_mas = pd.DataFrame.from_dict([ig.gen_function_tag(i, args.correo) for i in gen_id_mas.values()])
    dataframe_mas = pd.DataFrame.from_dict([res for res in (ig.gen_function_tag(i, args.correo) for i in gen_id_mas.values()) if res is not None])

    #dataframe_menos = pd.DataFrame.from_dict([ig.gen_function_tag(i, args.correo) for i in gen_id_menos.values()])
    dataframe_menos = pd.DataFrame.from_dict([res for res in (ig.gen_function_tag(i, args.correo) for i in gen_id_menos.values()) if res is not None])


    # Creacion de los archivos dataframe a csv
    
    logger.info('DataFrames listos')
    dataframe_mas.to_csv('../results/Genes_Sobreexpresados.csv', index=False)  
    dataframe_menos.to_csv('../results/Genes_Subexpresados.csv', index=False) 


