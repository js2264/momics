import configparser
import os
from pathlib import Path
from typing import Literal, Optional, Union
import tiledb


from .logging import logger

DEFAULT_CONFIG_PATH = Path.home() / ".momics.ini"
_has_logged = False


class LocalConfig:
    """
    LocalConfig is a class to manage the configuration of local momics repositories.
    """

    def __init__(self, path: Union[Path, str] = DEFAULT_CONFIG_PATH) -> None:
        """Initialize and read from a local config file."""

        if isinstance(path, str):
            path = Path(path)

        self.config_path = path
        self.cfg = configparser.ConfigParser()

        if self.config_path.exists():
            self.cfg.read(self.config_path)
        else:
            self.cfg.add_section("s3")
            self.cfg.add_section("gcs")
            self.cfg.add_section("azure")
            self._write_default_config()

    def _write_default_config(self) -> None:
        with open(self.config_path, "w") as configfile:
            self.cfg.write(configfile)

    def get(self, section, key, default=None) -> Union[str, None]:
        """Get value from the config file."""
        return self.cfg.get(section, key, fallback=default)

    def set(self, section, key, value) -> None:
        """Set value in the config file."""
        if not self.cfg.has_section(section):
            self.cfg.add_section(section)
        self.cfg.set(section, key, value)
        with open(self.config_path, "w") as configfile:
            self.cfg.write(configfile)

    def _validate_local_config(self) -> Literal[True]:
        """
        Validate LocalConfig by checking if all required keys exist for
        each cloud provider.
        """
        required_keys = {
            "s3": ["region", "access_key_id", "secret_access_key"],
            "gcs": ["project_id", "credentials"],
            "azure": ["account_name", "account_key"],
        }
        for section in self.cfg.sections():
            if not all([self.cfg.get(section, key, fallback=False) for key in required_keys[section]]):
                raise ValueError(f"Invalid cloud configuration. " f"Please provide all required values: {required_keys[section]}")
        return True


class S3Config:
    """
    This class is used to provide manual configuration for S3 access. It requires
    the following parameters:

    - `region`: The region of the S3 bucket
    - `access_key_id`: The access key ID for the S3 bucket
    - `secret_access_key`: The secret access key for the S3 bucket

    """

    def __init__(self, region=None, access_key_id=None, secret_access_key=None) -> None:
        self.region = region
        self.access_key_id = access_key_id
        self.secret_access_key = secret_access_key
        if not self._is_valid():
            raise ValueError(
                "Invalid S3 configuration. Please provide all required " "values: `region`, `access_key_id`, `secret_access_key`"
            )

    def _is_valid(self) -> bool:
        return all([self.region, self.access_key_id, self.secret_access_key])


class GCSConfig:
    """
    This class is used to provide manual configuration for GCS access. It requires
    the following parameters:

    - `project_id`: The project ID for the GCS bucket
    - `credentials`: Path to the credentials JSON file for the GCS bucket
    """

    def __init__(self, project_id=None, credentials=None) -> None:
        self.project_id = project_id
        self.credentials = credentials
        if not self._is_valid():
            raise ValueError("Invalid GCS configuration. Please provide all " "required values: `project_id`, `credentials`")

    def _is_valid(self) -> bool:
        return all([self.project_id, self.credentials])


class AzureConfig:
    """
    This class is used to provide manual configuration for Azure access. It requires
    the following parameters:

    - `account_name`: An account name
    - `account_key`: The associated account key
    """

    def __init__(self, account_name=None, account_key=None) -> None:
        self.account_name = account_name
        self.account_key = account_key
        if not self._is_valid():
            raise ValueError("Invalid Azure configuration. " "Please provide all required values: `account_name`, `account_key`")

    def _is_valid(self) -> bool:
        return all([self.account_name, self.account_key])


class MomicsConfig:
    """
    MomicsConfig is a class to automatically manage the configuration of the momics
    repository being accessed. It is used to set up the configuration for seamless
    cloud access.

    By default, it will parse a config file located at `~/.momics.ini` and
    use the first cloud configuration provided there. If the file does not exist or is
    incomplete, it will fall back to a base configuration which will only
    allow local momics repositories to be accessed.

    For example, a local config file might look like:

    ::

        [s3]
        region = us-west-1
        access_key_id = key
        secret_access_key = secret

        [gcs]
        region = us-west-1
        access_key_id = key
        secret_access_key = secret

        [azure]
        account_name = account
        account_key = key

    A manual configuration can also be used to access a specific cloud provider.

      - To access an S3 bucket, pass an `S3Config` object to the `s3` parameter.
      - To access a GCS bucket, pass a `GCSConfig` object to the `gcs` parameter.
      - To access an Azure bucket, pass an `AzureConfig` object to the `azure` \
          parameter.

    Cloud configuration is chosen in this order:
    `manual` (`S3Config` > `GCSConfig` > `AzureConfig`) > `local` > `blank`.
    """

    def __init__(
        self,
        s3: Optional[S3Config] = None,
        gcs: Optional[GCSConfig] = None,
        azure: Optional[AzureConfig] = None,
        local_cfg: Path = DEFAULT_CONFIG_PATH,
    ) -> None:
        """Initialize the momics configurator to enable cloud access."""

        # If a manual configuration is passed, use it.
        self.type: Optional[str] = None
        global _has_logged

        if s3 is not None:
            self.type = "s3"
            self.cfg = self._create_manual_tiledb_config(s3)
            if not _has_logged:
                _has_logged = True
                logger.info("Using S3 config.")
        elif gcs is not None:
            self.type = "gcs"
            self.cfg = self._create_manual_tiledb_config(gcs)
            if not _has_logged:
                _has_logged = True
                logger.info("Using GCS config.")
        elif azure is not None:
            self.type = "azure"
            self.cfg = self._create_manual_tiledb_config(azure)
            if not _has_logged:
                _has_logged = True
                logger.info("Using Azure config.")

        # Otherwise, parse local config.
        else:
            lcfg = LocalConfig(local_cfg)
            # If the file configuration is valid, use it.
            try:
                lcfg._validate_local_config()
                self.type = "local"
                self.cfg = self._create_tiledb_config_from_local(lcfg)
                logger.debug(f"Using local config from {lcfg.config_path} file.")
            # Otherwise, use blank configuration.
            except ValueError:
                self.cfg = tiledb.Config()
                if not _has_logged:
                    _has_logged = True
                    logger.info(
                        "No cloud config found for momics."
                        "Consider populating `~/.momics.ini` file with "
                        "configuration settings for cloud access."
                    )

        self.ctx = tiledb.cc.Context(self.cfg)
        self.vfs = tiledb.VFS(config=self.cfg, ctx=self.ctx)

    def _create_manual_tiledb_config(self, config: Union[S3Config, GCSConfig, AzureConfig]) -> tiledb.Config:
        if isinstance(config, S3Config):
            return tiledb.Config(
                {
                    "vfs.s3.region": config.region,
                    "vfs.s3.aws_access_key_id": config.access_key_id,
                    "vfs.s3.aws_secret_access_key": config.secret_access_key,
                }
            )
        elif isinstance(config, GCSConfig):
            os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = config.credentials
            return tiledb.Config(
                {
                    "vfs.gcs.project_id": config.project_id,
                    "vfs.gcs.credentials": config.credentials,
                }
            )
        elif isinstance(config, AzureConfig):
            return tiledb.Config(
                {
                    "vfs.azure.storage_account_name": config.account_name,
                    "vfs.azure.storage_account_key": config.account_key,
                }
            )

    def _create_tiledb_config_from_local(self, local_cfg) -> tiledb.Config:
        json = local_cfg.get("gcs", "credentials", default=None)
        if json is not None:
            os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = local_cfg.get("gcs", "credentials", default=None)
        return tiledb.Config(
            {
                "vfs.s3.region": local_cfg.get("s3", "region", default=None),
                "vfs.s3.aws_access_key_id": local_cfg.get("s3", "access_key_id", default=None),
                "vfs.s3.aws_secret_access_key": local_cfg.get("s3", "secret_access_key", default=None),
                #
                #
                "vfs.gcs.project_id": local_cfg.get("gcs", "project_id", default=None),
                "vfs.gcs.credentials": local_cfg.get("gcs", "credentials", default=None),
                #
                #
                "vfs.azure.storage_account_name": local_cfg.get("azure", "account_name", default=None),
                "vfs.azure.storage_account_key": local_cfg.get("azure", "account_key", default=None),
            }
        )
