#!/usr/bin/env python
#
# Initializes a git repository with some suggested best practices for DNAnexus
# WDL workflow development & continuous integration. Run
#    dxWDL_ci_init.py hello_world
# to create a subdirectory with that name, initialized as a local git repo,
# with the following which you can then customize:
#
# hello_world.wdl
#   Trivial WDL workflow template.
#
# hello_world_test.wdl
#   A testing workflow which runs the main workflow and inspects its results.
#
# build_workflow.py
#   A script that builds the workflow using dxWDL (downloading the dxWDL jar
#   file if needed), and optionally builds and runs the test workflow as well.
#
# DNAnexus platform project for continuous integration of the workflow.
# It has the following structure:
#   /builds/ - build_workflow.py by default builds the applets and workflow
#              into a subfolder named based on timestamp and git revision
#   /inputs/ - holds large and/or slow-evolving files or other data that may
#              be used as default inputs in the workflow and/or inputs to the
#              test.
# Data should never be deleted from the inputs folder, but moved into an 'old'
# subfolder when no longer actively used. This ensures that all historical
# versions of the workflow remain working, so long as they refer to these by
# unique ID. A git_revision property is set on the workflow.
#
# .travis.yml
# Travis CI script that installs dx and then runs build_workflow --test.
# A Travis CI secure environment variable with a long-lived DNAnexus API token
# must be added for this to work. Exact instructions will be shown when you run
# dxWDL_ci_init.py. The Travis CI integration is optional; you can always just
# build_workflow.py --test locally.

from __future__ import print_function
import dxpy
import argparse
import os
import sys
import errno
import subprocess

def main():
    argparser = argparse.ArgumentParser(description="Initialize a git repository for DNAnexus workflow development & continuous integration.")
    argparser.add_argument("name", help="workflow name (required)")
    argparser.add_argument("--project", help="DNAnexus project ID for continuous integration (default: create a new project)")
    args = argparser.parse_args()

    # initialize local git repository
    dir="./"+args.name
    if os.path.exists(dir):
        print("Local directory {} already exists. Aborting.".format(dir), file=sys.stderr)
        sys.exit(1)
    os.mkdir(dir)
    os.chdir(dir)
    subprocess.check_call(["git", "init"])

    # create the DNAnexus project for continuous integration
    if not args.project:
        args.project = dxpy.api.project_new({"name": "dxWDL_ci_{}".format(args.name)})["id"]
    project = dxpy.DXProject(args.project)
    initialize_project(project)

    # generate local applet and workflow source
    generate_wdl(args.name)
    generate_test_wdl(args.name)
    generate_test_input_json(args.name)
    generate_build_workflow_py(args.name, project)
    generate_travis_yml()
    generate_gitignore()

    # make initial git commit
    files_to_commit = [
        ".gitignore",
        args.name+".wdl",
        args.name+"_test.wdl",
        args.name+"_test.input.json",
        "build_workflow.py",
        ".travis.yml"
    ]
    subprocess.check_call(["git", "add"] + files_to_commit)
    subprocess.check_call(["git", "commit", "-m", "dxWDL_ci_init " + args.name])

    print("""
Initialized templates for dxWDL workflow continuous integration in {dir}
Modify them so '{dir}/build_workflow.py --test' does everything you'd like.
If you push the new repo to GitHub, you can hook it up to Travis CI with these
steps:
1. In Travis CI's Accounts settings, flip this repo's switch to ON.
2. Create a long-lived DNAnexus API token with CONTRIBUTE access to the
   continuous integration project (xxxx).
3. Run 'travis encrypt DX_AUTH_TOKEN=xxxx --add' in {dir}.
   ('gem install travis' first, if needed)
4. Commit the modified .travis.yml and push everything to GitHub.
Travis CI integration is optional; you can always just run
  {dir}/build_workflow.py --test
locally instead.
""".format(dir=dir))

def initialize_project(project):
    project.new_folder("/builds", parents=True)
    project.new_folder("/inputs/old", parents=True)
    print("Initialized DNAnexus project: {} ({})".format(project.describe()["name"], project.get_id()))

def generate_wdl(workflow_name):
    outfn=workflow_name+".wdl"
    with open(outfn, "w") as outf:
        print(hello_world_wdl_template.replace("<<NAME>>",workflow_name), file=outf)
    print("Generated example WDL in " + outfn)

hello_world_wdl_template = """workflow <<NAME>> {
  String who

  call hello_world {
    input:
      who = who
  }

  output {
    String message = hello_world.message
  }
}

task hello_world {
  String who

  command {
    set -ex -o pipefail
    echo "Hello, ${who}!" | tee message
  }

  output {
    String message = read_string("message")
  }
}
"""

def generate_test_wdl(workflow_name):
    outfn=workflow_name+"_test.wdl"
    with open(outfn, "w") as outf:
        print(test_wdl_template.replace("<<NAME>>",workflow_name), file=outf)
    print("Generated test WDL in " + outfn)

test_wdl_template = """import "<<NAME>>.wdl" as main

workflow <<NAME>>_test {  
  call main.<<NAME>>
  call assert_strings_equal {
    input:
      string1 = <<NAME>>.message
  }
}

task assert_strings_equal {
  String string1
  String string2

  command {
    set -ex -o pipefail
    if [ "${string1}" != "${string2}" ]; then
        echo "FAIL"
        exit 1
    fi
    echo "PASS"
  }
}
"""

def generate_test_input_json(workflow_name):
    outfn=workflow_name+"_test.input.json"
    with open(outfn, "w") as outf:
        print(test_input_json_template.replace("<<NAME>>",workflow_name), file=outf)
    print("Generated test input in " + outfn)
test_input_json_template = """{
    "<<NAME>>_test.<<NAME>>.who": "dxWDL",
    "<<NAME>>_test.assert_strings_equal.string2": "Hello, dxWDL!"
}
"""

def generate_build_workflow_py(workflow_name, project):
    with open("build_workflow.py", "w") as build_workflow_py:
        src = build_workflow_py_template.replace("<<NAME>>",workflow_name).replace("<<PROJECT>>",project.get_id())
        print(src, file=build_workflow_py)
    make_executable("build_workflow.py")
    print("Generated ./build_workflow.py")

build_workflow_py_template = """#!/usr/bin/env python
from __future__ import print_function
import dxpy
import argparse
import sys
import os
import subprocess
import json
import time

dxWDL_version="0.60.2"
here = os.getcwd()
git_revision = subprocess.check_output(["git", "describe", "--always", "--dirty", "--tags"], cwd=here).strip()

def main():
    argparser = argparse.ArgumentParser(description="Build <<NAME>> workflow on DNAnexus.")
    argparser.add_argument("--project", help="DNAnexus project ID", default="<<PROJECT>>")
    argparser.add_argument("--folder", help="Folder within project (default: timestamp/git-based)", default=None)
    argparser.add_argument("--test", help="Build and run test workflow", action='store_true')
    argparser.add_argument("--no-wait", help="With --test, launch the analysis and exit without awaiting completion", action='store_true')
    args = argparser.parse_args()

    ensure_dxWDL()

    # set up environment
    if args.folder is None:
        args.folder = time.strftime("/builds/%Y-%m-%d/%H%M%S-") + git_revision

    project = dxpy.DXProject(args.project)
    print("project: {} ({})".format(project.name, args.project))
    project.new_folder(args.folder, parents=True)
    print("folder: {}".format(args.folder))

    # build the workflow
    wf = dxWDL("<<NAME>>.wdl", project, args.folder)
    print("workflow: {} ({})".format(wf.name, wf.get_id()))

    # build and run the test, if desired
    if args.test:
        test_folder=args.folder+"/test"
        print("test folder: {}".format(test_folder))
        project.new_folder(test_folder)
        twf = dxWDL("<<NAME>>_test.wdl", project, test_folder, reorg=False, inputs="<<NAME>>_test.input.json")
        print("test workflow: {} ({})".format(twf.name, twf.get_id()))
        run_cmd=[
            "dx", "run", twf.get_id(),
            "--destination", "{}:{}".format(project.get_id(), test_folder),
            "--name", "<<NAME>> {} test".format(git_revision),
            "-f", "<<NAME>>_test.input.dx.json",
            "-y"
        ]
        if args.no_wait:
            subprocess.check_call(run_cmd)
        else:
            noise = subprocess.Popen(["/bin/bash", "-c", "while true; do sleep 60; date; done"])
            run_cmd = run_cmd + ["--wait"]
            try:
                subprocess.check_call(run_cmd)
                print("success")
            finally:
                noise.kill()

# download the dxWDL jar file if necessary
def ensure_dxWDL():
    dxWDL_fullname = "dxWDL-{}.jar".format(dxWDL_version)
    dxWDL_local_path = os.path.join(here, dxWDL_fullname)
    if not os.path.exists(dxWDL_local_path):
        # download the jar file
        download_cmd = [
            "wget", "-nv",
            "https://github.com/dnanexus-rnd/dxWDL/releases/download/{}/{}".format(dxWDL_version, dxWDL_fullname),
            "-O",
            dxWDL_local_path]
        print(" ".join(download_cmd))
        subprocess.check_call(download_cmd)
    return dxWDL_local_path

def dxWDL(filename, project, folder, reorg=True, inputs=None):
    dxWDL_path = ensure_dxWDL()
    cmd = ["java", "-jar", dxWDL_path, "compile",
           os.path.join(here, filename),
           "--project", project.get_id(),
           "--folder", folder]
    if inputs:
        cmd = cmd + ["--inputs", inputs]
    if reorg:
        cmd = cmd + ["--reorg"]

    buf = subprocess.check_output(cmd)
    wfid = buf.strip()
    wf = dxpy.DXWorkflow(wfid, project.get_id())
    wf.set_properties({"git_revision": git_revision})
    return wf

if __name__ == '__main__':
    main()
"""

def generate_travis_yml():
    with open(".travis.yml", "w") as outfile:
        print(travis_yml, file=outfile)
    # TODO: it would be nice to do 'travis encrypt DX_AUTH_TOKEN=xxxx'
    # automatically. However the repository has to be switched ON in Travis
    # before that command will work.
    print("Generated .travis.yml")

travis_yml = """language: python
python:
  - 2.7
# prevent travis from doing pip install requirements.txt
install: true
script:
# disable Travis default virtualenv
- deactivate
# deploy dx-toolkit
- wget https://wiki.dnanexus.com/images/files/dx-toolkit-current-ubuntu-12.04-amd64.tar.gz
- tar zxf dx-toolkit-current-ubuntu-12.04-amd64.tar.gz
- source dx-toolkit/environment
# execute workflow builder script
- python build_workflow.py --test
"""

def generate_gitignore():
    with open(".gitignore", "w") as outfile:
        print(gitignore_template, file=outfile)
    print("Generated .gitignore")
gitignore_template = """dxWDL-*.jar
*_test.input.dx.json
"""

def make_executable(path):
    mode = os.stat(path).st_mode
    mode |= (mode & 0444) >> 2    # copy R bits to X
    os.chmod(path, mode)

if __name__ == '__main__':
    main()
