import click
import os
import configparser

from . import cli

# Define the config file path
CONFIG_FILE = os.path.join(os.path.expanduser("~"), ".momics_config.ini")


def load_config():
    config = configparser.ConfigParser()
    if os.path.exists(CONFIG_FILE):
        config.read(CONFIG_FILE)
    return config


def save_config(config):
    with open(CONFIG_FILE, "w") as configfile:
        config.write(configfile)


@cli.group()
@click.pass_context
def config(ctx):
    """Manage S3 configuration"""
    pass


@config.command()
@click.argument("key")
@click.argument("value")
@click.pass_context
def set(ctx, key: str, value: str):
    """Set a configuration value"""
    config = load_config()
    if "s3" not in config:
        config["s3"] = {}
    config["s3"][key] = value
    save_config(config)
    click.echo(f"Set {key} to {value}")


@config.command()
@click.argument("key")
@click.pass_context
def get(ctx, key: str):
    """Get a configuration value"""
    config = load_config()
    if "s3" in config and key in config["s3"]:
        value = config["s3"][key]
        click.echo(f"{key}: {value}")
    else:
        click.echo(f"{key} not found")


@config.command()
@click.pass_context
def list(ctx):
    """List all configuration values"""
    config = load_config()
    if "s3" in config:
        for key, value in config["s3"].items():
            click.echo(f"{key}: {value}")
    else:
        click.echo("No S3 configuration found")
