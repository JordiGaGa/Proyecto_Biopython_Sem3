# ===========================================================================
# =                            imports
# ===========================================================================
import logging
# ===========================================================================
# =                            functions
# ===========================================================================

def getlogger(*args,**kwargs):

    logger = logging.getLogger()
    if not logger.hasHandlers():    
        logger.setLevel(logging.INFO)

        console_handler = logging.StreamHandler() 
        console_handler.setLevel(logging.INFO) 
        # Crear un formateador y añadirlo al manejador 
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s') 
        console_handler.setFormatter(formatter) 
        # Añadir el manejador al logger 
        logger.addHandler(console_handler)
    return logger

# ===========================================================================
# =                            main
# ===========================================================================
if __name__ == "__main__":
    pass