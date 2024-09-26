import logging


class MomicsFormatter(logging.Formatter):
    def format(self, record):
        return super().format(record)


log_format = "momics :: %(levelname)s :: %(asctime)s :: %(message)s"
formatter = MomicsFormatter(log_format)

handler = logging.StreamHandler()
handler.setFormatter(formatter)

logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.addHandler(handler)
