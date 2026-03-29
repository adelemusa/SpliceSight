import pyensembl
import os

_ensembl_release = None

def get_ensembl(release: int = 109, species: str = "human"):
    global _ensembl_release
    if _ensembl_release is None:
        _ensembl_release = pyensembl.EnsemblRelease(release, species=species)
        
        # Effettuerà il download in background se non già cachetato (~minuti al primo avvio)
        try:
            _ensembl_release.genes() 
        except ValueError:
            print("Indicizzazione PyEnsembl DB necessaria. Forzatura del download...")
            _ensembl_release.download()
            _ensembl_release.index()
            print("Indicizzazione PyEnsembl DB completata.")
    return _ensembl_release
