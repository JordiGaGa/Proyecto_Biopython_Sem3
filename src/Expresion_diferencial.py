# ===========================================================================
# =                            imports
# ===========================================================================
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from logger import getlogger
import sys
import os

# ===========================================================================
# =                            functions
# ===========================================================================

def analisis_diferencial(table:str, samples:dict, p_value:int=0.005) -> pd.DataFrame:
    '''
    Realiza un análisis diferencial en los datos proporcionados. 
    Args: 
        table (str): Ruta al archivo de la tabla de datos. 
        samples (dict[str, list]): Un diccionario donde las claves son nombres de muestras ("control" e "states") y los valores son listas de datos (list).
    Returns: 
        pd.Datafreme: Retorna un dataframe con los valores de expresion base y diferencial a "states"
    '''
    
    logger = getlogger('Analisis diferencial')
    # Eliminamos el output de la libreria pydeseq2
    
    original_stdout = sys.stdout
    original_stderr = sys.stderr


    sys.stderr = open(os.devnull, 'w')
    sys.stdout = open(os.devnull, 'w')

    logger.info('Iniciando analisis diferencial')
    try:    
        logger.info('Leyendo matriz de conteos')
        # Leer la matriz de conteo
        count_matrix = pd.read_csv(table,index_col=0).T
        if count_matrix.empty:
            logger.critical("EL dataframe cargado esta vacio")
            raise
        count_matrix = count_matrix.loc[[sample for group in samples.values() for sample in (group if isinstance(group,list) else [group])]]
        logger.info('Parseando matriz de conteos')
       # Reenombrar las filas en el orden deseado
        count_matrix = count_matrix.rename(index={sample : (group + sample[-1]) for group,cases in samples.items() for sample in (cases if isinstance(cases,list) else [cases])})
        count_matrix = count_matrix.round().astype(int)
    except Exception as e:
        logger.critical(f'Error creando la matriz de conteos: {e}')
        raise
    
    logger.info('Empezando analisis diferencial')
    
    try:
        #Agregar comentario
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

    except Exception as e:
        raise(f'Error al realizar el analisis diferencial: {e}')
    logger.info('Empezando analisis estadistico')
    
    try:
        ds = DeseqStats(dds,contrast=['conditions','states','control'])
        ds.summary()
    except Exception as e:
        raise(f'Error al realizar el analisis estadistico: {e}')
    #Regresamos los valores de expresion con un p-value ajustado menor a 0.05 
    
    sys.stdout = original_stdout
    sys.stderr = original_stderr
    
    return ds.results_df.query(f'padj<={p_value}')['log2FoldChange']

# ===========================================================================
# =                            Main
# ===========================================================================

if __name__ == "__main__":
    pass
