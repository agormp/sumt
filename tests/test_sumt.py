from types import SimpleNamespace

import phylotreelib as pt
import pytest
import sumt


@pytest.fixture
def annotated_sumtree():
    tree = pt.Tree.from_string("((A:0.123456789,B:2):3,C:4);")
    for node in tree.nodes:
        for attribute, value in {
            "height": 1.23456789, "height_sd": 0.15, "height_median": 1.2,
            "clade_cred": 0.9, "extra_meta": {"height_95%_CI": (0.5, 2.0), "unrequested": 99},
        }.items():
            tree.set_node_attribute(node, attribute, value)
    for children in tree.child_dict.values():
        for branch in children.values():
            branch.label = "0.9"
            branch.posterior = 0.95
            branch.bipartition_cred = 0.8
            branch.length_sd = 0.25
            branch.length_median = 2.5
            branch.extra_meta = {"length_95%_CI": (0.5, 4.0), "unrequested": 99}
    return tree


@pytest.fixture
def printing_args(tmp_path):
    return SimpleNamespace(
        treetype="hip", caheight=True, cladeheight=False, biplen=False, noblen=False,
        trackci=True, nometa=False, nolabel=False, outformat="nexus",
        outbase=tmp_path / "summary", nowarn=True, actively_rooted=False,
    )


@pytest.fixture
def silent_output():
    return SimpleNamespace(info=lambda *args, **kwargs: None)


class TestPrintSumtree:

    @pytest.mark.parametrize("treetype", ["con", "all", "mcc", "mbc", "hip", "mrhip"])
    @pytest.mark.parametrize("mode", ["biplen", "cladeheight", "caheight", "noblen"])
    @pytest.mark.parametrize("outformat", ["newick", "nexus"])
    @pytest.mark.parametrize("trackci", [False, True])
    def test_summary_annotations(self, annotated_sumtree, printing_args, silent_output,
                                 treetype, mode, outformat, trackci):
        args = printing_args
        args.treetype = treetype
        args.outformat = outformat
        args.trackci = trackci
        for name in ("biplen", "cladeheight", "caheight", "noblen"):
            setattr(args, name, name == mode)
        tree = annotated_sumtree
        plain = tree.newick()
        summary = SimpleNamespace(ci_labels=("95%_CI",))

        sumt.print_sumtree(tree, summary, args, silent_output)

        result = args.outbase.with_suffix("." + treetype).read_text()
        assert tree.rooted is (treetype in ("mcc", "hip", "mrhip"))
        if outformat == "nexus":
            declaration = "[&R] " if treetype in ("mcc", "hip", "mrhip") else "[&U] "
            assert "tree nexus_tree = " + declaration in result
        else:
            assert "[&R]" not in result and "[&U]" not in result
        assert "posterior=" not in result
        assert "label=" not in result
        assert "unrequested=" not in result
        if treetype in ("mcc", "hip", "mrhip"):
            assert "clade_cred=0.9" in result
            assert ")0.9" in result
        else:
            assert ("bipartition_cred=0.8" in result) == (mode != "noblen")
            assert "clade_cred=" not in result
            assert ")0.8" in result
        if mode in ("cladeheight", "caheight"):
            assert "height=1.234568" in result
            assert "height_sd=0.15" in result
            assert ("height_median=1.2" in result) == trackci
            if trackci:
                assert "height_95%_CI={0.5,2.0}" in result
        elif mode == "biplen":
            assert "length_sd=0.25" in result
            assert ("length_median=2.5" in result) == trackci
            if trackci:
                assert "length_95%_CI={0.5,4.0}" in result
        else:
            assert ":" not in result
        if not trackci:
            assert "_CI=" not in result
        assert result.startswith("#NEXUS") == (outformat == "nexus")
        assert result.endswith("\n")
        assert tree.newick() == plain

    @pytest.mark.parametrize("outformat", ["newick", "nexus"])
    def test_metadata_and_label_suppression(self, annotated_sumtree, printing_args,
                                           silent_output, outformat):
        args = printing_args
        args.nometa = True
        args.nolabel = True
        args.outformat = outformat
        sumt.print_sumtree(annotated_sumtree, SimpleNamespace(ci_labels=("95%_CI",)),
                          args, silent_output)
        result = args.outbase.with_suffix(".hip").read_text()
        if outformat == "nexus":
            assert "tree nexus_tree = [&R] " in result
        assert "[&" not in result.replace("[&R]", "")
        assert ")0.9" not in result
        assert ":0.1234568" in result

    @pytest.mark.parametrize("treetype", ["con", "all", "mbc"])
    @pytest.mark.parametrize("actively_rooted", [False, True, ["A"], []])
    @pytest.mark.parametrize("nometa", [False, True])
    def test_explicit_rooting_declaration(self, annotated_sumtree, printing_args, silent_output,
                                         treetype, actively_rooted, nometa):
        args = printing_args
        args.treetype = treetype
        args.actively_rooted = actively_rooted
        args.nometa = nometa
        sumt.print_sumtree(annotated_sumtree, SimpleNamespace(ci_labels=("95%_CI",)),
                          args, silent_output)
        result = args.outbase.with_suffix("." + treetype).read_text()
        declaration = "[&R]" if actively_rooted else "[&U]"
        assert "tree nexus_tree = " + declaration + " " in result
        if nometa:
            assert "[&" not in result.replace(declaration, "")

    def test_no_interval_labels(self, annotated_sumtree, printing_args, silent_output):
        sumt.print_sumtree(annotated_sumtree, SimpleNamespace(ci_labels=()),
                          printing_args, silent_output)
        result = printing_args.outbase.with_suffix(".hip").read_text()
        assert "height=" in result
        assert "_CI=" not in result
        assert "height_median" not in result
