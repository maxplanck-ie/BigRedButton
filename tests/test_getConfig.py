import configparser

import pytest

from BRB import getConfig


def _writeConfig(path, extraLines=""):
    path.write_text(f"[Paths]\nbaseData=/tmp/base\n{extraLines}")


class TestGetConfigArgumentValidation:
    def test_missing_configfile_argument_exits(self, capsys):
        with pytest.raises(SystemExit) as excinfo:
            getConfig.getConfig(None)
        assert excinfo.value.code == 1
        assert "configFile is undefined" in capsys.readouterr().out

    def test_nonexistent_configfile_exits(self, tmp_path, capsys):
        missing = tmp_path / "does_not_exist.ini"
        with pytest.raises(SystemExit) as excinfo:
            getConfig.getConfig(str(missing))
        assert excinfo.value.code == 1
        assert "does not exists" in capsys.readouterr().out

    def test_missing_paths_section_exits(self, tmp_path, capsys):
        configFile = tmp_path / "brb.ini"
        configFile.write_text("[Options]\nsleepTime=1\n")
        with pytest.raises(SystemExit) as excinfo:
            getConfig.getConfig(str(configFile))
        assert excinfo.value.code == 1
        assert "No Paths defined" in capsys.readouterr().out


class TestGetConfigSuccess:
    def test_returns_configparser_with_paths_section(self, tmp_path, monkeypatch):
        configFile = tmp_path / "brb.ini"
        _writeConfig(configFile)
        monkeypatch.setattr(getConfig, "configGitInfo", lambda *a, **k: "deadbeef")

        config = getConfig.getConfig(str(configFile))

        assert isinstance(config, configparser.ConfigParser)
        assert config.get("Paths", "baseData") == "/tmp/base"

    def test_sets_options_configcommit_from_configgitinfo(self, tmp_path, monkeypatch):
        configFile = tmp_path / "brb.ini"
        _writeConfig(configFile)
        monkeypatch.setattr(getConfig, "configGitInfo", lambda *a, **k: "abc1234")

        config = getConfig.getConfig(str(configFile))

        assert config.get("Options", "configCommit") == "abc1234"

    def test_configcommit_falls_back_to_empty_string_when_no_git_info(
        self, tmp_path, monkeypatch
    ):
        configFile = tmp_path / "brb.ini"
        _writeConfig(configFile)
        monkeypatch.setattr(getConfig, "configGitInfo", lambda *a, **k: None)

        config = getConfig.getConfig(str(configFile))

        assert config.get("Options", "configCommit") == ""

    def test_creates_options_section_if_absent(self, tmp_path, monkeypatch):
        configFile = tmp_path / "brb.ini"
        _writeConfig(configFile)
        monkeypatch.setattr(getConfig, "configGitInfo", lambda *a, **k: "abc1234")

        config = getConfig.getConfig(str(configFile))

        assert "Options" in config.sections()

    def test_preserves_existing_options_section_entries(self, tmp_path, monkeypatch):
        configFile = tmp_path / "brb.ini"
        _writeConfig(configFile, extraLines="[Options]\nsleepTime=2\n")
        monkeypatch.setattr(getConfig, "configGitInfo", lambda *a, **k: "abc1234")

        config = getConfig.getConfig(str(configFile))

        assert config.get("Options", "sleepTime") == "2"
        assert config.get("Options", "configCommit") == "abc1234"

    def test_passes_custom_git_binary_from_software_section(
        self, tmp_path, monkeypatch
    ):
        configFile = tmp_path / "brb.ini"
        _writeConfig(configFile, extraLines="[software]\ngit=/usr/local/bin/git\n")
        seen = {}

        def fakeConfigGitInfo(configfile, gitBin):
            seen["gitBin"] = gitBin

        monkeypatch.setattr(getConfig, "configGitInfo", fakeConfigGitInfo)

        getConfig.getConfig(str(configFile))

        assert seen["gitBin"] == "/usr/local/bin/git"

    def test_defaults_git_binary_to_git_when_software_section_absent(
        self, tmp_path, monkeypatch
    ):
        configFile = tmp_path / "brb.ini"
        _writeConfig(configFile)
        seen = {}

        def fakeConfigGitInfo(configfile, gitBin):
            seen["gitBin"] = gitBin

        monkeypatch.setattr(getConfig, "configGitInfo", fakeConfigGitInfo)

        getConfig.getConfig(str(configFile))

        assert seen["gitBin"] == "git"
