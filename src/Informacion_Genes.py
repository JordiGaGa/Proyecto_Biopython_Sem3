# ===========================================================================
# =                            imports
# ===========================================================================

from Bio import Entrez
import time


# ===========================================================================
# =                            functions
# ===========================================================================

def Id_Gene(gen,org,co):
    """
    Funcion para obtener los IDs de los genes de una lista proporcionada
    
    Args:
        gen (string) : Value de la lista de genes que viene siendo el nombre del gen

        org (string) : Nombre del organismo a analizar

        co (string)  : Correo del usuario para acceden a Entrez
    Returns:
        record["IdList"] : Regresa el ID del gen en caso de encontrar, y si no regresa None

    """
    Entrez.email = co
    query = f"{gen}[Gene] AND {org}[Organism]"
    with Entrez.esearch(db="gene", term=query) as handle:
        record = Entrez.read(handle)
        time.sleep(.12)
    handle.close()  
    return record["IdList"][0] if record["IdList"] else None

def gen_function_tag(id,co):
    """
    Funcion para regresar dataframe que contenga:
    ID del gen, nombre del gen, descripción de la proteína y locus tag
    
    Args:
        id (string) : Value de una lista el cual es el ID de un dado gen

        co (string)  : Correo del usuario para acceden a Entrez
    Returns:
        info(pd.Dataframe) : DataFrame con los datos mencionados
        None : En caso de error se regresa un None 
    """
    Entrez.email = co
    # Busca y regresa la informacion del gen en formato xml 
    with Entrez.efetch(db="gene", id=id, retmode="xml") as handle:
        record = Entrez.read(handle)[0]
        time.sleep(.12)
        # Analizamos si el id del gen tiene todos los campos buscados
        try: 
            gn = record['Entrezgene_gene']['Gene-ref'].get('Gene-ref_locus', None),
            pt = record['Entrezgene_prot']['Prot-ref'].get('Prot-ref_desc', None),
            lt = record['Entrezgene_gene']['Gene-ref'].get('Gene-ref_locus-tag', None)
        # Si no cumple los campos se regresa un error 
        except Exception as err:
            handle.close()
            return
        
        # Se crea el diccionario con los datos extraidos
        info = {
            "gene_id": id,
            "gene_name":record['Entrezgene_gene']['Gene-ref'].get('Gene-ref_locus', None),
            "prot_desc":record['Entrezgene_prot']['Prot-ref'].get('Prot-ref_desc', None),
            "locus_tag":record['Entrezgene_gene']['Gene-ref'].get('Gene-ref_locus-tag', None)
        }

    handle.close()
    return(info)
