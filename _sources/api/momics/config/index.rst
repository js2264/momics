momics.config
=============

.. py:module:: momics.config


Attributes
----------

.. autoapisummary::

   momics.config.DEFAULT_CONFIG_PATH


Classes
-------

.. autoapisummary::

   momics.config.AzureConfig
   momics.config.GCSConfig
   momics.config.LocalConfig
   momics.config.MomicsConfig
   momics.config.S3Config


Module Contents
---------------

.. py:class:: AzureConfig(account_name=None, account_key=None)

   This class is used to provide manual configuration for Azure access. It requires
   the following parameters:

   - `account_name`: An account name
   - `account_key`: The associated account key


   .. py:attribute:: account_key
      :value: None



   .. py:attribute:: account_name
      :value: None



.. py:class:: GCSConfig(project_id=None, credentials=None)

   This class is used to provide manual configuration for GCS access. It requires
   the following parameters:

   - `project_id`: The project ID for the GCS bucket
   - `credentials`: Path to the credentials JSON file for the GCS bucket


   .. py:attribute:: credentials
      :value: None



   .. py:attribute:: project_id
      :value: None



.. py:class:: LocalConfig(path = DEFAULT_CONFIG_PATH)

   LocalConfig is a class to manage the configuration of local momics repositories.


   .. py:method:: get(section, key, default=None)

      Get value from the config file.



   .. py:method:: set(section, key, value)

      Set value in the config file.



   .. py:attribute:: cfg


   .. py:attribute:: config_path


.. py:class:: MomicsConfig(s3 = None, gcs = None, azure = None, local_cfg = DEFAULT_CONFIG_PATH)

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
     - To access an Azure bucket, pass an `AzureConfig` object to the `azure`           parameter.

   Cloud configuration is chosen in this order:
   `manual` (`S3Config` > `GCSConfig` > `AzureConfig`) > `local` > `blank`.


   .. py:attribute:: ctx


   .. py:attribute:: type
      :type:  Optional[str]
      :value: None



   .. py:attribute:: vfs


.. py:class:: S3Config(region=None, access_key_id=None, secret_access_key=None)

   This class is used to provide manual configuration for S3 access. It requires
   the following parameters:

   - `region`: The region of the S3 bucket
   - `access_key_id`: The access key ID for the S3 bucket
   - `secret_access_key`: The secret access key for the S3 bucket



   .. py:attribute:: access_key_id
      :value: None



   .. py:attribute:: region
      :value: None



   .. py:attribute:: secret_access_key
      :value: None



.. py:data:: DEFAULT_CONFIG_PATH

