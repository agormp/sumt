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
    tree.rooted = True
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
    return SimpleNamespace(info=lambda *args, **kwargs: None,
                           force=lambda *args, **kwargs: None)


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
        tree.rooted = treetype in ("mcc", "hip", "mrhip")
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
    @pytest.mark.parametrize("rooted", [None, False, True])
    @pytest.mark.parametrize("nometa", [False, True])
    def test_stored_rooting_declaration(self, annotated_sumtree, printing_args, silent_output,
                                         treetype, rooted, nometa):
        args = printing_args
        args.treetype = treetype
        args.actively_rooted = not rooted
        annotated_sumtree.rooted = rooted
        args.nometa = nometa
        sumt.print_sumtree(annotated_sumtree, SimpleNamespace(ci_labels=("95%_CI",)),
                          args, silent_output)
        result = args.outbase.with_suffix("." + treetype).read_text()
        declaration = "[&R]" if rooted is True else "[&U]" if rooted is False else ""
        assert annotated_sumtree.rooted is rooted
        assert "tree nexus_tree = " + (declaration + " " if declaration else "") in result
        if rooted is None:
            assert "[&R]" not in result and "[&U]" not in result
        if nometa:
            assert "[&" not in result.replace(declaration, "")

    def test_no_interval_labels(self, annotated_sumtree, printing_args, silent_output):
        sumt.print_sumtree(annotated_sumtree, SimpleNamespace(ci_labels=()),
                          printing_args, silent_output)
        result = printing_args.outbase.with_suffix(".hip").read_text()
        assert "height=" in result
        assert "_CI=" not in result
        assert "height_median" not in result


class TestComputeSumtree:

    @pytest.mark.parametrize("treetype", ["con", "all", "mcc", "mbc", "hip", "mrhip"])
    @pytest.mark.parametrize("actively_rooted", [False, True, ["A"], []])
    @pytest.mark.parametrize("cpus", [1, 2])
    def test_workflow_establishes_rootedness(self, monkeypatch, printing_args, silent_output,
                                            treetype, actively_rooted, cpus):
        args = printing_args
        args.treetype = treetype
        args.actively_rooted = actively_rooted
        args.cpus = cpus
        args.outgroup = args.rootmid = args.rootminvar = False
        tree = SimpleNamespace(rooted=None, set_blens_from_heights=lambda: None)
        monkeypatch.setattr(pt, "build_sumtree", lambda *a, **kw: tree)
        monkeypatch.setattr(sumt, "set_ca_heights_concurrent",
                            lambda *a, **kw: (tree, {42}))
        monkeypatch.setattr(pt, "TreePostProcessor",
                            lambda summary: SimpleNamespace(annotate_sumtree=lambda t: t))
        worker_pids = set()
        result = sumt.compute_sumtree(object(), args, [], silent_output, 1, worker_pids)
        assert result is tree
        assert result.rooted is (treetype in ("mcc", "hip", "mrhip") or bool(actively_rooted))
        assert worker_pids == ({42} if cpus == 2 else set())


class TestFastTreecount:
    def test_lexical_comments_and_root_suffix(self, tmp_path):
        path = tmp_path / "annotated.nexus"
        path.write_text(
            '#NEXUS\nbegin trees;\n'
            'TREE a=[&R](A[&note="tree fake = (;);"]:1,B:2)ROOT:0.5;\n'
            'UTREE b=(A:1,B:2)[&height=3]:0.0;\nend;\n')
        args = SimpleNamespace(informat="nexus")
        assert sumt.fast_treecount(path, args) == 2


class TestChunkedTreeStringsFromFiles:
    def test_context_and_translation_blocks_stay_with_descriptions(self, tmp_path):
        path = tmp_path / "mixed.nexus"
        path.write_text(
            '#NEXUS\nbegin trees;\ntranslate 1 A, 2 B;\n'
            'TREE first=[&R](1[&note="a;b"]:1,2:2)ROOT:0.5;\n'
            'UTREE second=(1:1,2:2)SECOND:0.25;\nend;\n'
            'begin trees;\ntranslate 1 C, 2 D;\n'
            'TREE third=(1:1,2:2)THIRD:0.0;\nend;\n')
        chunks = list(sumt.chunked_tree_strings_from_files([(3, 0, path)], 100))
        assert [len(chunk) for chunk, _, _ in chunks] == [2, 1]
        first_chunk, file_idx, first_parser = chunks[0]
        second_chunk, _, second_parser = chunks[1]
        assert file_idx == 0 and first_parser is not second_parser
        assert '[&note="a;b"]' in first_chunk[0][0]
        assert first_chunk[0][0].endswith('ROOT:0.5;')
        assert [kind for _, kind in first_chunk] == ["TREE", "UTREE"]
        assert [kind for _, kind in second_chunk] == ["TREE"]
        first = [pt.Tree._from_string_private(first_parser, text, command_kind=kind)
                 for text, kind in first_chunk]
        second = [pt.Tree._from_string_private(second_parser, text, command_kind=kind)
                  for text, kind in second_chunk]
        assert [tree.rooted for tree in first + second] == [True, False, None]
        assert all(tree.leaves == {"A", "B"} for tree in first)
        assert second[0].leaves == {"C", "D"}

    def test_burnin_skips_without_losing_next_context(self, tmp_path):
        path = tmp_path / "burnin.nexus"
        path.write_text('#NEXUS\nbegin trees;\nTREE a=[&R](A:1,B:2);\n'
                        'UTREE b=(A:1,B:2)ROOT:0.5;\nend;\n')
        [(chunk, _, parser)] = list(sumt.chunked_tree_strings_from_files([(2, 1, path)], 10))
        assert len(chunk) == 1 and chunk[0][1] == "UTREE"
        tree = pt.Tree._from_string_private(parser, chunk[0][0], command_kind=chunk[0][1])
        assert tree.rooted is False and tree.root_suffix.label == "ROOT"
        assert tree.root_suffix.length == 0.5


def _parser_worker_context_fixture(tmp_path):
    path = tmp_path / "workers.nexus"
    path.write_text('#NEXUS\nbegin trees;\nTREE a=[&R](A[&note="x;y"]:1,B:2);\n'
                    'UTREE b=(A:1,B:2);\nTREE c=(A:1,B:2);\nend;\n')
    [(chunk, _, parser)] = list(sumt.chunked_tree_strings_from_files([(3, 0, path)], 10))
    return chunk, parser


class TestWorkerProcessChunk:
    def test_intact_text_and_per_description_rootedness(self, tmp_path, monkeypatch):
        chunk, parser = _parser_worker_context_fixture(tmp_path)
        seen = []
        original = pt.Tree._from_string_private.__func__
        def record(cls, parser_obj, text, **kwargs):
            tree = original(cls, parser_obj, text, **kwargs)
            seen.append(tree.rooted)
            return tree
        monkeypatch.setattr(pt.Tree, "_from_string_private", classmethod(record))
        def forbidden(*args, **kwargs):
            raise AssertionError("worker stripped complete annotated input")
        monkeypatch.setattr(pt, "remove_comments", forbidden)
        args = SimpleNamespace(trackbips=True, trackclades=True, trackblen=False,
            trackroot=False, trackrootblen=False, trackheight=False, tracktopo=False,
            track_subcladepairs=False, treeprobs=False, trackci=False, ci_probs=(), ci_k=7)
        summary, file_idx, pid = sumt.worker_process_chunk(chunk, 42, parser, args)
        assert len(summary) == 3 and file_idx == 42 and pid > 0
        assert seen == [True, False, None]


class TestWorkerProcessCaChunk:
    def test_intact_text_and_per_description_rootedness(self, tmp_path, monkeypatch):
        chunk, parser = _parser_worker_context_fixture(tmp_path)
        plan = pt.CAHeightEstimator.build_plan(pt.Tree.from_string("(A:1,B:2);"))
        monkeypatch.setattr(sumt, "_CA_PLAN", plan)
        monkeypatch.setattr(sumt, "_CA_TRACKCI", False)
        monkeypatch.setattr(sumt, "_CA_QUANTILE_K", 7)
        seen = []
        original = pt.Tree._from_string_private.__func__
        def record(cls, parser_obj, text, **kwargs):
            tree = original(cls, parser_obj, text, **kwargs)
            seen.append(tree.rooted)
            return tree
        monkeypatch.setattr(pt.Tree, "_from_string_private", classmethod(record))
        def forbidden(*args, **kwargs):
            raise AssertionError("CA worker stripped complete annotated input")
        monkeypatch.setattr(pt, "remove_comments", forbidden)
        estimator, count, pid = sumt.worker_process_ca_chunk(chunk, parser)
        assert count == 3 and pid > 0
        assert all(record.n == 3 for record in estimator.acc.values())
        assert seen == [True, False, None]


@pytest.fixture
def topology_probability_summary():
    trees = [pt.Tree.from_string(text, parse_meta=True) for text in (
        '((A[&height=0]:1,B:1)AB[&height=1]:1,(C:1,D:1):1)ROOT:0.5[&length_sd=0.1];',
        '((A:1,C:1):1,(B:1,D:1):1);',
        '((A:1,D:1):1,(B:1,C:1):1);',
    )]
    transdict = trees[0].transdict()
    return SimpleNamespace(
        trackclades=True, trackbips=False, transdict=transdict,
        translateblock=trees[0].translateblock(transdict),
        cladetoposummary={i: SimpleNamespace(posterior=p, tree=tree)
                         for i, (p, tree) in enumerate(zip((0.25, 0.5, 0.25), trees))},
        biptoposummary={i: SimpleNamespace(posterior=p, tree=tree)
                       for i, (p, tree) in enumerate(zip((0.5, 0.25, 0.25), trees))},
    )


class TestComputeTrprobs:

    @pytest.mark.parametrize('trackclades,trackbips', [(True, False), (False, True), (True, True)])
    def test_distribution_precedence(self, topology_probability_summary, trackclades, trackbips):
        summary = topology_probability_summary
        summary.trackclades, summary.trackbips = trackclades, trackbips
        args = SimpleNamespace(trackclades=not trackclades, trackbips=not trackbips)
        selected = summary.cladetoposummary if trackclades else summary.biptoposummary
        expected = sorted(selected.values(), key=lambda entry: entry.posterior, reverse=True)

        result = sumt.compute_trprobs(summary, args)

        assert [p for p, tree in result] == [0.5, 0.25, 0.25]
        assert all(tree is entry.tree for (_, tree), entry in zip(result, expected))


class TestPrintTrprobs:

    @pytest.mark.parametrize('trackclades,trackbips,rooted,actively_rooted,treetype,declaration', [
        (True, False, False, False, 'mcc', '[&R]'),
        (False, True, True, True, 'con', '[&U]'),
        (True, True, None, False, 'con', '[&R]'),
    ])
    def test_distribution_declaration_preserves_representatives(
            self, topology_probability_summary, printing_args, silent_output, monkeypatch,
            trackclades, trackbips, rooted, actively_rooted, treetype, declaration):
        summary, args = topology_probability_summary, printing_args
        summary.trackclades, summary.trackbips = trackclades, trackbips
        args.trackclades, args.trackbips = not trackclades, not trackbips
        args.actively_rooted, args.treetype, args.treeprobs = actively_rooted, treetype, 1.0
        entries = sumt.compute_trprobs(summary, args)
        before = []
        for _, tree in entries:
            tree.rooted = rooted
            before.append((tree.root, tree.root_suffix, tree.newick(),
                           tuple((n, id(owner)) for n, owner in tree.nodedict.items()),
                           tuple((p, c, id(branch)) for p, children in tree.child_dict.items()
                                 for c, branch in children.items())))

        def forbidden(*args, **kwargs):
            pytest.fail('Printing topology probabilities must not copy or reroot representatives')

        for method in ('copy_treeobject', 'reroot', 'unroot'):
            monkeypatch.setattr(pt.Tree, method, forbidden)
        sumt.print_trprobs(summary, entries, args, silent_output)

        path = args.outbase.with_suffix('.trprobs')
        result = path.read_text()
        lines = [line for line in result.splitlines() if line.startswith('    tree ')]
        for i, ((freq, tree), state, line) in enumerate(zip(entries, before, lines), 1):
            assert line == (f'    tree tree_{i} [p = {freq:.6f}] '
                            f'[P = {sum(p for p, _ in entries[:i]):.6f}] = {declaration} '
                            + tree.newick(printdist=False, printlabels=False,
                                          transdict=summary.transdict))
            assert tree.rooted is rooted
            assert (tree.root, tree.root_suffix, tree.newick(),
                    tuple((n, id(owner)) for n, owner in tree.nodedict.items()),
                    tuple((p, c, id(branch)) for p, children in tree.child_dict.items()
                          for c, branch in children.items())) == state
        assert len(lines) == len(entries)
        assert summary.translateblock in result
        with pt.Nexustreefile(path) as reader:
            parsed = list(reader)
        assert len(parsed) == len(entries)
        assert all(tree.rooted is trackclades for tree in parsed)
        assert [tree.topology() for tree in parsed] == [tree.topology() for _, tree in entries]

    @pytest.mark.parametrize('cutoff,expected_count', [(0.6, 2), (0.75, 3), (1.0, 3)])
    def test_probability_headers_and_strict_cumulative_cutoff(
            self, topology_probability_summary, printing_args, silent_output, cutoff, expected_count):
        summary, args = topology_probability_summary, printing_args
        args.treeprobs = cutoff
        entries = sumt.compute_trprobs(summary, args)

        sumt.print_trprobs(summary, entries, args, silent_output)

        result = args.outbase.with_suffix('.trprobs').read_text()
        if cutoff < 1:
            header = (f'[This file contains the {round(cutoff * 100)}% most probable trees found during the\n'
                      f'MCMC search, sorted by posterior probability (the {round(cutoff * 100)}% HPD interval).\n')
        else:
            header = ('[This file contains all trees that were found during the MCMC\n'
                      'search, sorted by posterior probability. \n')
        expected = ('#NEXUS\n\n' + header
                    + "Lower case 'p' indicates the posterior probability of a tree.\n"
                    + "Upper case 'P' indicates the cumulative posterior probability.]\n\n"
                    + 'begin trees;\n' + summary.translateblock)
        for i, cumulative in enumerate((0.5, 0.75, 1.0)[:expected_count], 1):
            freq, tree = entries[i - 1]
            expected += (f'    tree tree_{i} [p = {freq:.6f}] [P = {cumulative:.6f}] = '
                         + tree.newick(printdist=False, printlabels=False,
                                       transdict=summary.transdict) + '\n')
        expected += 'end;\n'
        assert result.replace('= [&R] ', '= ') == expected
