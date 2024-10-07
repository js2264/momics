import configparser
import os

import click

from . import cli

# Define the config file path
CONFIG_FILE = os.path.join(os.path.expanduser("~"), ".momics.ini")


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
    """Manage cloud provider configurations"""
    pass


def cloud_provider_group(name):
    @config.group(name=name)
    @click.pass_context
    def cloud_provider(ctx):
        """Manage configuration for cloud provider."""
        pass

    @cloud_provider.command()
    @click.argument("key")
    @click.argument("value")
    @click.pass_context
    def set(ctx, key: str, value: str):
        """Set a configuration value"""
        config = load_config()
        if name not in config:
            config[name] = {}
        config[name][key] = value
        save_config(config)
        click.echo(f"Set {key} to {value} for {name}")

    @cloud_provider.command()
    @click.argument("key")
    @click.pass_context
    def get(ctx, key: str):
        """Get a configuration value"""
        config = load_config()
        if name in config and key in config[name]:
            value = config[name][key]
            click.echo(f"{key}: {value} for {name}")
        else:
            click.echo(f"{key} not found in {name} configuration")

    @cloud_provider.command()
    @click.pass_context
    def list(ctx):
        """List all configuration values"""
        config = load_config()
        if name in config:
            for key, value in config[name].items():
                click.echo(f"{key}: {value} for {name}")
        else:
            click.echo(f"No configuration found for {name}")


# Add subcommands for each cloud provider
cloud_provider_group("s3")
cloud_provider_group("gcs")
cloud_provider_group("azure")
