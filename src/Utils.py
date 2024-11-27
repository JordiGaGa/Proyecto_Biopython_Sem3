import pandas as pd
import numpy as np

def nivel_exp_df(df: pd.DataFrame, conteo: bool = False):
  """
  Filtra y agrega una nueva columna de acuerdo con los cambios en el nivel de
  expresion.

  Args:
        df (pd.DataFrame): El dataframe que contiene los datos con una columna llamada
                           'log2FoldChange' que representa los cambios en
                           los niveles de expresión.
  Returns:
        pd.Dataframe: El dataframe ingresado con una nueva columna que contenga el cambio
                      en el nivel de expresion de cada gen, ya sea que este sobreexpresado
                      subexpresado o que no haya cambios significativos dependiendo del umbral
                      ingresado.
  """
  # Sacar copia del df para no modificar el original
  exp_level = df.copy()

  # Clasificar los niveles de expresión y agregar una nueva columna 'expresion_change'
  exp_level['expresion_change'] = exp_level['log2FoldChange'].apply(
                                  lambda x: 'Subexpresado' if x < 0 else 'Sobreexpresado')

  # Devolver el conteo dependiente del cambio de expresion
  if conteo:
    return exp_level['expresion_change'].value_counts().to_dict()

  # Devolver el dataframe modificado
  return exp_level


def outliers_dif_exp(df_exp: pd.DataFrame):
  """
  Encuentra los valores de los outliers de la distribucion de genes sobre y sub expresados

  Args:
        df_exp (pd.DataFrame): El dataframe original que contiene la información de expresión de genes.

    Returns:
        dict: Diccionario cuyas llaves son el nivel de expresion ('Subexpresado' o 'Sobreexpresado') y los
              valores corresponden a una lista que contenga cada uno de los outliers de cada distribucion
  """
  df_sobreexpresados = df_exp[df_exp['log2FoldChange'] > 0]
  df_subexpresados = df_exp[df_exp['log2FoldChange'] < 0]

  for data in [df_sobreexpresados,df_subexpresados]:
    Q1 = np.percentile(data['log2FoldChange'], 25)
    Q3 = np.percentile(data['log2FoldChange'], 75)
    IQR = Q3 - Q1
    if data['log2FoldChange'].iloc[0] > 0:
      limite = Q3 + 1.5 * IQR
      sobreexp_outliers = list(data[data['log2FoldChange'] > limite].index)
    else:
      limite = Q1 - 1.5 * IQR
      subexp_outliers = list(data[data['log2FoldChange'] < limite].index)


  return {'Sobreexpresado':sobreexp_outliers,'Subexpresado':subexp_outliers }