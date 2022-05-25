from hmm2parquet import ProfilesFile

def test_parser():
    tests = [
        "testfiles/single_profile.hmm"
    ]
    for test in tests:
        file = ProfilesFile()
        file.read_file(test)
        file.to_parquet("testfiles")