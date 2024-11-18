import pandas as pd
from pydeseq2.dds import DeseqDataSet
import argparse 
import os

def process_list(argline):
    if os.path.isfile(argline):
        with open(argline) as file:
            return [ int(txid) for line in file.readlines() if (txid:=line.strip())]
    else:
        return [str(name) for name in argline.split(',') if name]

def analisis_diferencial(table:str, samples:dict) -> pd.DataFrame:
    '''
    Realiza un análisis diferencial en los datos proporcionados. 
    Args: 
        table (str): Ruta al archivo de la tabla de datos. 
        samples (dict[str, list]): Un diccionario donde las claves son nombres de muestras ("control" e "states") y los valores son listas de datos (list).
    Returns: 
        pd.Datafreme: Retorna un dataframe con los valores de expresion base y diferencial a "states"
    '''
    try:    
        # Leer la matriz de conteo
        count_matrix = pd.read_csv(table,index_col=0).T
        if count_matrix.empty:
            raise ValueError("EL dataframe cargado esta vacio")
        count_matrix = count_matrix.loc[[sample for group in samples.values() for sample in (group if isinstance(group,list) else [group])]]

        # Reenombrar las filas en el orden deseado
        count_matrix = count_matrix.rename(index={sample : (group + sample[-1]) for group,cases in samples.items() for sample in (cases if isinstance(cases,list) else [cases])})
        count_matrix = count_matrix.round().astype(int)
    except Exception as e:
        raise(f'Error creando la matriz de conteos: {e}')
    try:
        # Definir la metadata con "low-mg" como la condición de referencia
        metadata = pd.DataFrame({
            'conditions' : [case[:-1] for case in list(count_matrix.index)] 
        }, index=list(count_matrix.index))
        
        # Crear el objeto DeseqDataSet
        dds = DeseqDataSet(
            counts=count_matrix,
            metadata=metadata,
            design_factors="conditions"
        )

        # Ejecutar el análisis DESeq2
        dds.deseq2()

        # Tabla Log natural
        dds.varm['LFC']
    except Exception as e:
        raise(f'Error al realizar el analisis diferenical: {e}')
    return dds.varm['LFC']

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Recibe una matrix de conteos y regresa un dataframe con el analisis de expresion')
    
    parser.add_argument('-i','--input', type=str, required=True,
                        help='Ruta a la matriz de conteo')
    parser.add_argument('-cn','--control',type=process_list, required=True,
                        help='Lista con el nombre de las columnas que son el control en el dataframe')
    parser.add_argument('-st','--states',type=process_list, required=True,
                        help='Lista con el nombre de las columnas que son el estado alternativo en el dataframe')
    parser.add_argument('-sv','--save', action='store_true',
                        help='Salvar o no el dataframe final')
    
    args = parser.parse_args()

    samples = {'states':args.states, 'control':args.control}
    print((df:=analisis_diferencial(table=args.input,samples=samples)))
    if args.save and not df.empty:
        df.to_csv('analisis.txt')
