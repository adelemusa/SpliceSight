import pytest
import os
import sys
import tempfile
from io import StringIO

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "backend"))

from core.rmats_parser import (
    EVENT_TYPES,
    detect_event_type,
    get_file_type,
    parse_inc_level,
    check_proteome_impact,
    extract_coordinates,
    parse_single_file,
    parse_rmats_files,
    RMatsConfig,
    SplicingEvent,
    get_available_event_types,
    _count_by_type,
)


class TestDetectEventType:
    def test_se_detection(self):
        assert detect_event_type("SE.MATS.JC.txt") == "SE"
        assert detect_event_type("SE.MATS.JCEC.txt") == "SE"
        assert detect_event_type("fromGTF.SE.txt") == "SE"

    def test_a3ss_detection(self):
        assert detect_event_type("A3SS.MATS.JC.txt") == "A3SS"
        assert detect_event_type("A3SS.MATS.JCEC.txt") == "A3SS"

    def test_a5ss_detection(self):
        assert detect_event_type("A5SS.MATS.JC.txt") == "A5SS"
        assert detect_event_type("A5SS.MATS.JCEC.txt") == "A5SS"

    def test_mxe_detection(self):
        assert detect_event_type("MXE.MATS.JC.txt") == "MXE"
        assert detect_event_type("MXE.MATS.JCEC.txt") == "MXE"

    def test_ri_detection(self):
        assert detect_event_type("RI.MATS.JC.txt") == "RI"
        assert detect_event_type("RI.MATS.JCEC.txt") == "RI"

    def test_unknown_type(self):
        assert detect_event_type("summary.txt") is None
        assert detect_event_type("random.txt") is None
        assert detect_event_type("") is None


class TestGetFileType:
    def test_jc_file(self):
        assert get_file_type("SE.MATS.JC.txt") == "JC"

    def test_jcec_file(self):
        assert get_file_type("SE.MATS.JCEC.txt") == "JCEC"


class TestParseIncLevel:
    def test_valid_values(self):
        assert parse_inc_level("0.5,0.6,0.7") == [0.5, 0.6, 0.7]
        assert parse_inc_level("1.0,0.5,NA,0.3") == [1.0, 0.5, None, 0.3]

    def test_empty_and_na(self):
        assert parse_inc_level("") == []
        assert parse_inc_level("NA") == []
        assert parse_inc_level(None) == []


class TestCheckProteomeImpact:
    def test_frameshift(self):
        impact, is_nmd = check_proteome_impact("TEST", 100)
        assert "Frameshift" in impact
        assert is_nmd is True

    def test_in_frame(self):
        impact, is_nmd = check_proteome_impact("TEST", 99)
        assert "In-frame" in impact
        assert is_nmd is False

    def test_terminal_exon(self):
        impact, is_nmd = check_proteome_impact("TEST", 100, is_terminal=True)
        assert "Frameshift" in impact
        assert is_nmd is False


class TestExtractCoordinates:
    def test_se_coordinates(self):
        row = {
            "exonStart_0base": 1000,
            "exonEnd": 1100,
            "upstreamES": 900,
            "upstreamEE": 950,
            "downstreamES": 1150,
            "downstreamEE": 1200,
        }
        cols = EVENT_TYPES["SE"]["exon_cols"]
        coords = extract_coordinates(row, "SE", cols)
        assert coords["start"] == 1000
        assert coords["end"] == 1100

    def test_missing_coordinates(self):
        row = {"exonStart_0base": 1000}
        cols = EVENT_TYPES["SE"]["exon_cols"]
        coords = extract_coordinates(row, "SE", cols)
        assert "start" in coords
        assert "end" not in coords


class TestCountByType:
    def test_count_events(self):
        events = [
            SplicingEvent("1", "SE", "G1", "ENSG1", "chr1", "+", 0.01, 0.1, 0.05),
            SplicingEvent("2", "SE", "G2", "ENSG2", "chr1", "+", 0.02, 0.2, 0.04),
            SplicingEvent("3", "A3SS", "G3", "ENSG3", "chr1", "+", 0.01, 0.3, 0.03),
            SplicingEvent("4", "RI", "G4", "ENSG4", "chr1", "+", 0.01, 0.4, 0.02),
        ]
        counts = _count_by_type(events)
        assert counts == {"SE": 2, "A3SS": 1, "RI": 1}

    def test_empty_events(self):
        counts = _count_by_type([])
        assert counts == {}


class TestRMatsConfig:
    def test_default_config(self):
        config = RMatsConfig()
        assert config.fdr_threshold == 0.05
        assert config.dpsi_threshold == 0.1
        assert config.include_jcec is True
        assert config.species == "human"
        assert config.ensembl_release == 109

    def test_custom_config(self):
        config = RMatsConfig(fdr_threshold=0.01, dpsi_threshold=0.2)
        assert config.fdr_threshold == 0.01
        assert config.dpsi_threshold == 0.2


class TestGetAvailableEventTypes:
    def test_returns_all_types(self):
        types = get_available_event_types()
        type_names = [t["type"] for t in types]
        assert "SE" in type_names
        assert "A3SS" in type_names
        assert "A5SS" in type_names
        assert "MXE" in type_names
        assert "RI" in type_names

    def test_has_descriptions(self):
        types = get_available_event_types()
        for t in types:
            assert "description" in t
            assert "patterns" in t


class TestParseSingleFile:
    @pytest.fixture
    def se_test_file(self):
        header = "ID\tGeneID\tgeneSymbol\tchr\tstrand\texonStart_0base\texonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\tID\tIJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\tIncFormLen\tSkipFormLen\tPValue\tFDR\tIncLevel1\tIncLevel2\tIncLevelDifference"
        row1 = "1\tENSG000001\tgeneA\tchr1\t+\t1000\t1100\t900\t950\t1150\t1200\t1\t10,10\t5,5\t5,5\t5,5\t98\t49\t0.001\t0.01\t0.5,0.5\t0.6,0.6\t-0.1"
        row2 = "2\tENSG000002\tgeneB\tchr2\t-\t2000\t2200\t1900\t1950\t2250\t2300\t2\t20,20\t10,10\t10,10\t10,10\t98\t49\t0.001\t0.03\t0.7,0.7\t0.5,0.5\t0.2"
        row3 = "3\tENSG000003\tgeneC\tchr3\t+\t3000\t3100\t2900\t2950\t3150\t3200\t3\t5,5\t2,2\t2,2\t2,2\t98\t49\t0.001\t0.1\t0.4,0.4\t0.5,0.5\t-0.1"
        content = header + "\n" + row1 + "\n" + row2 + "\n" + row3 + "\n"

        tmpdir = tempfile.gettempdir()
        file_path = os.path.join(tmpdir, "SE.MATS.JC.txt")
        with open(file_path, "w") as f:
            f.write(content)
        yield file_path
        if os.path.exists(file_path):
            os.unlink(file_path)

    def test_parse_se_file(self, se_test_file):
        config = RMatsConfig(fdr_threshold=0.05, dpsi_threshold=0.1)
        events, genes = parse_single_file(se_test_file, config)
        assert len(events) == 1
        assert "geneB" in genes

    def test_parse_empty_file(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("ID\tFDR\tIncLevelDifference\n")
            f.flush()
            path = f.name

        try:
            config = RMatsConfig()
            events, genes = parse_single_file(path, config)
            assert len(events) == 0
            assert len(genes) == 0
        finally:
            os.unlink(path)


class TestParseRmatsFiles:
    @pytest.fixture
    def test_directory(self):
        content_se = """ID	GeneID	geneSymbol	chr	strand	exonStart_0base	exonEnd	upstreamES	upstreamEE	downstreamES	downstreamEE	ID	IJC_SAMPLE_1	SJC_SAMPLE_1	IJC_SAMPLE_2	SJC_SAMPLE_2	IncFormLen	SkipFormLen	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference
1	ENSG001	gene1	chr1	+	1000	1100	900	950	1150	1200	1	10	5	5	98	49	0.001	0.01	0.5	0.6	-0.1
2	ENSG002	gene2	chr2	+	2000	2200	1900	1950	2250	2300	2	20	10	10	98	49	0.001	0.03	0.7	0.5	0.2
"""
        content_a3ss = """ID	GeneID	geneSymbol	chr	strand	longExonStart_0base	longExonEnd	shortES	shortEE	flankingES	flankingEE	ID	IJC_SAMPLE_1	SJC_SAMPLE_1	IJC_SAMPLE_2	SJC_SAMPLE_2	IncFormLen	SkipFormLen	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference
1	ENSG003	gene3	chr3	+	1000	1100	1000	1050	1150	1200	1	10	5	5	72	49	0.001	0.02	0.5	0.6	-0.1
"""
        with tempfile.TemporaryDirectory() as tmpdir:
            se_path = os.path.join(tmpdir, "SE.MATS.JC.txt")
            a3ss_path = os.path.join(tmpdir, "A3SS.MATS.JC.txt")

            with open(se_path, "w") as f:
                f.write(content_se)
            with open(a3ss_path, "w") as f:
                f.write(content_a3ss)

            yield tmpdir

    def test_parse_directory(self, test_directory):
        config = RMatsConfig(fdr_threshold=0.05, dpsi_threshold=0.1)
        result = parse_rmats_files(test_directory, config=config)

        assert "events" in result
        assert "summary" in result
        assert result["summary"]["total_events"] >= 0

    def test_parse_nonexistent_path(self):
        result = parse_rmats_files("/nonexistent/path")
        assert result["events"] == []
        assert result["summary"]["total_events"] == 0


class TestSplicingEvent:
    def test_event_creation(self):
        event = SplicingEvent(
            event_id="1",
            event_type="SE",
            gene_symbol="BRCA1",
            gene_id="ENSG00001",
            chromosome="17",
            strand="+",
            fdr=0.01,
            dpsi=0.5,
            p_value=0.001,
            coordinates={"start": 1000, "end": 1100},
            proteome_impact="Frameshift (100bp)",
            nmd_candidate=True,
        )
        assert event.event_id == "1"
        assert event.event_type == "SE"
        assert event.nmd_candidate is True


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
