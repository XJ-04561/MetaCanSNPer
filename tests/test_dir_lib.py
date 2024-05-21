
from MetaCanSNPer.core.DirectoryLibrary import DirectoryLibrary

def test_filename_alignment():
    """Testing alignment scenarios."""
    from MetaCanSNPer.core.DirectoryLibrary import PathList
    
    aligned1 = PathList("FSC458_R1.fq", "FSC458_R2.fq")

    assert "FSC458" == aligned1.name, '"FSC458" == PathList("FSC458_R1.fq", "FSC458_R2.fq").name'

    aligned2 = PathList("FSC458.fq", "FSC458(1).fq")

    assert "FSC458" == aligned2.name, '"FSC458" == PathList("FSC458.fq", "FSC458(1).fq").name'
    
    aligned3 = PathList("FSC458.fna", "FSC658-[FSC458_R1].fq", "FSC567-[FSC458_R2].fq")

    assert "FSC458" == aligned3.name, 'PathList("FSC458.fna", "FSC658-[FSC458_R1].fq", "FSC567-[FSC458_R2].fq")'

def test_appdir_mirroring():
    
    from appdirs import AppDirs

    AD = AppDirs("MetaCanSNPer", multipath=True)
    DL = DirectoryLibrary("francisella_tularensis", ["FSC458_R1.fq.gz", "FSC458_R2.fq.gz"])

    assert DL.userDataDir == AD.user_data_dir
    assert DL.siteDataDir == AD.site_data_dir
    assert DL.userConfigDir == AD.user_config_dir
    assert DL.siteConfigDir == AD.site_config_dir
    assert DL.userCacheDir == AD.user_cache_dir
    assert DL.userLogDir == AD.user_log_dir

def test_initialization():
    dl1 = DirectoryLibrary("francisella_tularensis", "FSC458.fq")
    dl2 = DirectoryLibrary("francisella_tularensis", ["FSC458_R1.fq", "FSC458_R2.fq"])

    assert dl1.queryName == dl2.queryName