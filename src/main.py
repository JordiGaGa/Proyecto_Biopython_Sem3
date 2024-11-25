# ===========================================================================
# =                            imports
# ===========================================================================
from Expresion_diferencial import analisis_diferencial
import argparse
import os
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
                        help='Ruta a la matriz de conteo')
    parser.add_argument('-cn','--control',type=process_list, required=True,
                        help='Lista con el nombre de las columnas que son el control en el dataframe')
    parser.add_argument('-st','--states',type=process_list, required=True,
                        help='Lista con el nombre de las columnas que son el estado alternativo en el dataframe')
    parser.add_argument('-sv_all','--save_full_data', action='store_true',
                        help='Salvar o no el dataframe final')
    
    args, _ = parser.parse_known_args()

    if args.save_full_data:
        parser.add_argument('-o_all','--output_full_data',default='./full_data.csv',
                            help='Directorio donde guardar el archivo con todos los datos de expresion')
        args = parser.parse_args()

    try:
        samples = {'states':args.states, 'control':args.control}
        
        df = analisis_diferencial(table=args.input,samples=samples)
        #Añadir Resto de Procesos
       
        if args.save_full_data and not df.empty:
            df.to_csv(args.output_full_data)
        
    except Exception as e:
        raise(f'Putamadre{e}')
