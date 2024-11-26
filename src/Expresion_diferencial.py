# ===========================================================================
# =                            imports
# ===========================================================================
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from logger import getlogger
import argparse 
import sys
import os

# ===========================================================================
# =                            functions
# ===========================================================================

def analisis_diferencial(table_csv:str, samples:dict) -> pd.DataFrame:
    '''
    Realiza un análisis diferencial en los datos proporcionados.

    Args: 
        table_csv (str): Ruta al archivo de la tabla de datos en formato csv. 
        samples (dict[str, list]): Un diccionario donde las llaves son el tipo de muestras, ya sea 
                                   que estas pertenescan a un control ("control") o que sean las distintas
                                   condiciones ("states") y los valores son la lista de nombres de las muestras
                                   de las respectivas muestras.
    Returns: 
        pd.Datafreme: Retorna un dataframe con los valores de expresion base y diferencial a "states"
    '''
    
    logger = getlogger('Analisis diferencial')

    # Guardar referencia a los outputs
    original_stdout = sys.stdout
    original_stderr = sys.stderr

    # Suprimir el output impreso a pantalla (por parte de las librerias)
    sys.stderr = open(os.devnull, 'w')
    sys.stdout = open(os.devnull, 'w')

    logger.info('Iniciando analisis diferencial')
    try:    
        logger.info('Leyendo matriz de conteos')
        # Leer matriz de conteo
        count_matrix = pd.read_csv(table_csv,index_col=0).T
        if count_matrix.empty:
            logger.critical("EL dataframe cargado esta vacio")
            raise
        count_matrix = count_matrix.loc[[sample for group in samples.values() for sample in (group if isinstance(group,list) else [group])]]
        logger.info('Parseando matriz de conteos')
       # Reenombrar filas en el orden deseado
        count_matrix = count_matrix.rename(index={sample : (group + sample[-1]) for group,cases in samples.items() for sample in (cases if isinstance(cases,list) else [cases])})
        count_matrix = count_matrix.round().astype(int)
    except Exception as e:
        logger.critical(f'Error creando la matriz de conteos: {e}')
        raise
    
    logger.info('Empezando analisis deferencial')
    
    try:
        # Generar archivo de metadata para el proceso de analisis diferencial por pydefse2
        metadata = pd.DataFrame({
            'conditions' : [case[:-1] for case in list(count_matrix.index)] 
        }, index=list(count_matrix.index))
        # Crear objeto DeseqDataSet
        dds = DeseqDataSet(
            counts=count_matrix,
            metadata=metadata,
            design_factors="conditions"
        )
        # Ejecutar análisis DESeq2
        dds.deseq2()
    except Exception as e:
        raise(f'Error al realizar el analisis diferenical: {e}')
    

    logger.info('Empezando analisis estadistico')
    
    try:
        ds = DeseqStats(dds,contrast=['conditions','states','control'])
        ds.summary()
    except Exception as e:
        raise(f'Error al realizar el analisis estadistico: {e}')
    
    # Activar nuevamente el output
    sys.stdout = original_stdout
    sys.stderr = original_stderr
    
    # Regresar valores de expresion con un p-value ajustado menor a 0.05 
    return ds.results_df.query('padj<=0.05')['log2FoldChange']

# ===========================================================================
# =                            Main
# ===========================================================================

if __name__ == "__main__":
    pass
