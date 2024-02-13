import os

TEST_ENV = os.environ.get("TEST_ENV", default="local")

if TEST_ENV == "local":
    import sys

    sys.path.append("../duet_tools")
