import pandas as pd
from pydeseq2.dds import DeseqDataSet

def analisis_diferencial(tabla:str = ""):
    # Leer la matriz de conteo
    count_matrix = pd.read_csv(path,index_col=0).T
    count_matrix = count_matrix.drop('length')

    # Reorganizar las filas en el orden deseado
    count_matrix = count_matrix.reindex(['ns1', 'ns2', 'low-mg1', 'low-mg2'])
    count_matrix = count_matrix.rename(index={'ns1': 'control1', 'ns2': 'control2'})
    count_matrix = count_matrix.round().astype(int)

    # Definir la metadata con "low-mg" como la condición de referencia
    metadata = pd.DataFrame({
        'condition': ['control', 'control', 'low-mg', 'low-mg']  # Invertimos la posición de "low-mg" y "ns"
    }, index=['control1', 'control2', 'low-mg1', 'low-mg2'])

    # Crear el objeto DeseqDataSet
    dds = DeseqDataSet(
        counts=count_matrix,
        metadata=metadata,
        design_factors="condition"
    )

    # Ejecutar el análisis DESeq2
    dds.deseq2()

    # Tabla Log natural
    dds.varm['LFC']

    return dds.varm['LFC']

if __name__ == "__main__":
    return