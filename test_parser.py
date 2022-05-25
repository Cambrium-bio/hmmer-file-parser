from hmm2parquet import ProfilesFile

def test_parser():
    tests = [
        "Pfam-A.hmm"
    ]
    for test in tests:
        file = ProfilesFile()
        file.read_file(test)
        file.to_parquet("pfam_output")