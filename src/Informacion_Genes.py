# ===========================================================================
# =                            imports
# ===========================================================================

from Bio import Entrez
import pandas as pd
import argparse
import time
from Expresion_diferencial import analisis_diferencial, process_list

# ===========================================================================
# =                            functions
# ===========================================================================
def Lista_genes(dds,num_genes):
    """
    Funcion para regresar 2 listas, una con los genes más sobreexpresados y otra
    con los más subexpresados basado en la cantidad de genes que indique el usuario
    
    Args:
        dds (pd.DataFrame): El dataframe que contiene los datos con una columna llamada 
                           'conditions_states_vs_control' que representa los cambios en 
                           los niveles de expresión.
        num_genes (int): Cantidad de genes que el usuario quiere que regrese la función

    Returns:
        genes_mas (list) : Lista con los genes más sobreexpresados
        genes menos (list) : Lista con los genes más subexpresados

    """
    # Creamos los sets de los genes mas sobreexpresados y mas subexpresados 
    top10 = dds['conditions_states_vs_control'].nlargest(num_genes)
    gen_mas = set(top10.index)
    bottom10 = dds['conditions_states_vs_control'].nsmallest(num_genes)
    gen_menos = set(bottom10.index)

    # Convertimos los sets a listas 
    genes_mas = [*gen_mas]
    genes_menos = [*gen_menos]

    return genes_mas, genes_menos

def Id_Gene(gen,org):
    """
    Funcion para obtener los IDs de los genes de una lista proporcionada
    
    Args:
        gen (string) : Value de la lista de genes que viene siendo el nombre del gen

        org (string) : Nombre del organismo a analizar
    Returns:
        record["IdList"] : Regresa el ID del gen en caso de encontrar, y si no regresa None

    """
    query = f"{gen}[Gene] AND {org}[Organism]"
    with Entrez.esearch(db="gene", term=query) as handle:
        record = Entrez.read(handle)
        time.sleep(.5)
    handle.close()  
    return record["IdList"][0] if record["IdList"] else None

def gen_function_tag(id):
    """
    Funcion para regresar dataframe que contenga:
    ID del gen, nombre del gen, descripción de la proteína y locus tag
    
    Args:
        id (string) : Value de una lista el cual es el ID de un dado gen
    Returns:
        info(pd.Dataframe) : DataFrame con los datos mencionados
    """
    # Busca y regresa la informacion del gen en formato xml 
    with Entrez.efetch(db="gene", id=id, retmode="xml") as handle:
        record = Entrez.read(handle)[0]
        time.sleep(.5)
    # Accede al campo en la estructura de datos
    info = {
        "gene_id": id,
        "gene_name":record['Entrezgene_gene']['Gene-ref'].get('Gene-ref_locus', None),
        "prot_desc":record['Entrezgene_prot']['Prot-ref'].get('Prot-ref_desc', None),
        "locus_tag":record['Entrezgene_gene']['Gene-ref'].get('Gene-ref_locus-tag', None)
        #### checar si hay mas info
    }
    handle.close()
    return(info)

# ===========================================================================
# =                            Main
# ===========================================================================

# Variables dadas por el usuario
parser = argparse.ArgumentParser(description='Regresa un 2 dataframes de los genes sobreexpresados y subexpresados'
                                'con su ID nombre del gen, descripcion de la proteina y locus tag')
    
parser.add_argument('-co', '--correo', type=str, required=True,
                        help='Correo electronico para acceder a Entrez') # jordigg@lcg.unam.mx 
parser.add_argument('-i','--input', type=str, required=True, 
                        help='Ruta a la matriz de conteo') # '../data/GSE276379_RNASeq_kallisto.csv'
parser.add_argument('-cn','--control',type=process_list, required=True,
                        help='Lista con el nombre de las columnas que son el control en el dataframe') # n1,n2
parser.add_argument('-st','--states',type=process_list, required=True,
                        help='Lista con el nombre de las columnas que son el estado alternativo en el dataframe') # low-mg1,low-mg2
parser.add_argument('-ng', '--numgenes', type=int, required=True, 
                        help='Cantidad de genes sobreexpresados y subexpresados que deseas que se regresen') # 10
parser.add_argument('-org', '--organismo', type=str, required=True,
                        help='Nombre del organismo del cual se hizo el analisis de expresion diferencial') # 'escherichia coli'
 
args = parser.parse_args()

Entrez.email = args.correo 

samples = {'states':args.states, 'control':args.control}

#archivo = "../data/GSE276379_RNASeq_kallisto.csv"
#samples = {'states':['low-mg1', 'low-mg2'], 'control':['ns1', 'ns2']}
#Entrez.email = "aggonzal@lcg.unam.mx"
#num_genes = 10 
#organismo = "escherichia coli"


# Llamada a función de análisis diferencial 

dds = analisis_diferencial(table=args.input,samples= samples)

# Llamada a función que regresa las listas de genes sobreexpresados y subexpresados
genes_mas,genes_menos = Lista_genes(dds, args.numgenes)

# Llamada a función que regresa los IDs de los genes a buscar
gen_id_mas = {gen:Id_Gene(gen, args.organismo) for gen in genes_mas}
gen_id_menos = {gen:Id_Gene(gen, args.organismo) for gen in genes_menos}

# Llamada a función que crea los dataframes con ID, gen, descripción de la proteína y locus tag
dataframe_mas = pd.DataFrame.from_dict([gen_function_tag(i) for i in gen_id_mas.values()])
dataframe_menos = pd.DataFrame.from_dict([gen_function_tag(i) for i in gen_id_menos.values()])

# Resultado impreso en pantalla
print(dataframe_mas)
print(dataframe_menos)