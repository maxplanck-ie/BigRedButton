import configparser

import pytest

from BRB import getConfig


def _writeConfig(path, extraLines=""):
    path.write_text(
        f"[Paths]\nbaseData=/tmp/base\n[Options]\nvalidAnalysisTypes=dna\n{extraLines}"
    )


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
        configFile.write_text("[Options]\nvalidAnalysisTypes=dna\nsleepTime=1\n")
        with pytest.raises(SystemExit) as excinfo:
            getConfig.getConfig(str(configFile))
        assert excinfo.value.code == 1
        assert "No Paths defined" in capsys.readouterr().out

    def test_legacy_options_key_exits(self, tmp_path, capsys):
        configFile = tmp_path / "brb.ini"
        configFile.write_text(
            "[Paths]\nbaseData=/tmp/base\n[Options]\nvalidLibraryTypes=dna\n"
        )

        with pytest.raises(SystemExit) as excinfo:
            getConfig.getConfig(str(configFile))

        assert excinfo.value.code == 1
        output = capsys.readouterr().out
        assert "validLibraryTypes" in output
        assert "validAnalysisTypes" in output
        assert str(configFile) in output

    def test_legacy_external_key_exits(self, tmp_path, capsys):
        configFile = tmp_path / "brb.ini"
        configFile.write_text(
            "[Paths]\nbaseData=/tmp/base\n[Options]\nvalidAnalysisTypes=dna\n[external]\nLibraryTypes=dna\n"
        )

        with pytest.raises(SystemExit) as excinfo:
            getConfig.getConfig(str(configFile))

        assert excinfo.value.code == 1
        output = capsys.readouterr().out
        assert "[external] LibraryTypes" in output
        assert "AnalysisTypes" in output

    def test_missing_valid_analysis_types_exits(self, tmp_path, capsys):
        configFile = tmp_path / "brb.ini"
        configFile.write_text("[Paths]\nbaseData=/tmp/base\n[Options]\nsleepTime=1\n")

        with pytest.raises(SystemExit) as excinfo:
            getConfig.getConfig(str(configFile))

        assert excinfo.value.code == 1
        assert "validAnalysisTypes" in capsys.readouterr().out


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

    def test_valid_analysis_types_config_loads(self, tmp_path, monkeypatch):
        configFile = tmp_path / "brb.ini"
        _writeConfig(configFile)
        monkeypatch.setattr(getConfig, "configGitInfo", lambda *a, **k: "abc1234")

        config = getConfig.getConfig(str(configFile))

        assert "Options" in config.sections()
        assert config.get("Options", "validAnalysisTypes") == "dna"

    def test_preserves_existing_options_section_entries(self, tmp_path, monkeypatch):
        configFile = tmp_path / "brb.ini"
        _writeConfig(configFile, extraLines="sleepTime=2\n")
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
