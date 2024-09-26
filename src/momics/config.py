from pathlib import Path
import configparser
import tiledb


class MomicsConfig:
    """
    MomicsConfig is a class that manages the configuration of the momics
    repository being accessed. It is used to set up the configuration for seamless
    cloud access.
    """

    def __init__(self, path=Path.home() / ".momics_config.ini"):
        """Initialize the momics configurator to cloud access."""

        self.config_path = path
        if not self.config_path.exists():
            self.config_path.touch()

        cfg = configparser.ConfigParser()
        cfg.read(self.config_path)

        if len(cfg.sections()) > 0:
            if "s3" in cfg.sections():
                tiledb_cfg = tiledb.Config(
                    {
                        "vfs.s3.region": cfg["s3"]["region"],
                        "vfs.s3.aws_access_key_id": cfg["s3"]["access_key_id"],
                        "vfs.s3.aws_secret_access_key": cfg["s3"]["secret_access_key"],
                    }
                )
            else:
                tiledb_cfg = tiledb.Config()
        else:
            tiledb_cfg = tiledb.Config()

        self.cfg = tiledb_cfg
        self.ctx = tiledb.cc.Context(tiledb_cfg)
        self.vfs = tiledb.VFS(config=self.cfg, ctx=self.ctx)
