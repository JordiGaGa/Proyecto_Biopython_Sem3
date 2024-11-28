import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from Utils import nivel_exp_df

def pie_expresion_plot(df: pd.Series, save : bool=False, output_dir: str = '.'):
    """
    Genera un gráfico de pie con los cambios en los niveles de expresion mostrados en el
    dataframe ingresado. Si se introduce un directorio de salida, la imagen sera
    guardada dentro en formato '.png'.

    Args:
        df (pd.DataFrame): El dataframe que contiene los datos con una columna llamada
                           'conditions_states_vs_control' que representa los cambios en
                           los niveles de expresión.
        output_dir (str, opcional): Ruta pcon el directorio para guardar la imagen
                                    generada en formato '.png'.


    Returns:
        plt.Figure: La figura del gráfico de pie generada.
    """

    # Clasificar los valores en categorías
    categorias = nivel_exp_df(df,True)

    # Crear gráfico de pie
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.pie(categorias.values(),
           wedgeprops={"linewidth": 2, "edgecolor": "white"},
           autopct='%1.1f%%',
           startangle=140)
    ax.set_title('Distribucion del cambio de expresion en genes con un p-valor significativo')
    ax.axis('equal')
    ax.legend(categorias.keys(),
              loc="center left",
              bbox_to_anchor=(1, 0, 0.5, 1))

    fig.tight_layout()

    # Guardar la imagen si se proporciona una ruta
    if save:
        plt.savefig(f'{output_dir}/pie_expresion_plot.jpg', transparent=True, bbox_inches= 'tight')
        print(f'Imagen "pie_expresion_plot.png" guardada dentro del directorio {output_dir}')

    return fig

def expresion_dist_plot(df: pd.Series, graph: str = 'box', save : bool = False, output_dir: str = '.'):
  """
  Genera un gráfico de puntos con los cambios en los niveles de expresion mostrados en el
  dataframe ingresado. Si se introduce un directorio de salida, la imagen sera
  guardada dentro en formato '.png'.

  Args:
      df (pd.DataFrame): El dataframe que contiene los datos con una columna llamada
                          'log2FoldChange' que representa los cambios en
                          los niveles de expresión.
      output_dir (str, opcional): Ruta pcon el directorio para guardar la imagen
                                  generada en formato '.png'.


  Returns:
      plt.Figure: Box plot con los valores de log2FoldChange de los genes sobre y sub expresados.
  """
  # Copiar el dataframe para evitar modificar el original
  exp_level = nivel_exp_df(df)

  # Crear figura y agregar un grid
  f, ax = plt.subplots(figsize=(7, 6))
  ax.xaxis.grid(True, linestyle='--', linewidth=0.5)

  # Crear boxplot o violinplot
  if graph == 'box':
    sns.boxplot(data=exp_level, x='log2FoldChange', y='expresion_change', hue ='expresion_change', ax=ax)
  elif graph == 'vl':
    sns.violinplot(data=exp_level, x='log2FoldChange', y='expresion_change', hue ='expresion_change', ax=ax)

  ax.set(ylabel=None)  # Eliminar el label del eje Y
  ax.set_title('Valores de cambio con p-value significativo')

  f.tight_layout()

  # Guardarlo si es necesario
  if save:
        plt.savefig(f'{output_dir}/expresion_plot.jpg', transparent=True)
        print(f'Imagen "expresion_plot.png" guardada dentro del directorio {output_dir}')

  # Mostrar la imagen
  return f