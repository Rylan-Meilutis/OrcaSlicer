"""Exercise the publisher's actual shell gate without network or artifacts."""

import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


WORKFLOW = Path(__file__).resolve().parents[2] / ".github/workflows/publish_release.yml"


@unittest.skipUnless(shutil.which("bash") and shutil.which("jq"), "bash and jq required")
class ReleaseBuildGateTest(unittest.TestCase):
    def check_gate(self, *, accepted=False, annotated=False, tag_sha="abc",
                   wait=False, pending=False, **overrides):
        source = WORKFLOW.read_text()
        step = source.split("      - name: Verify source build run\n", 1)[1]
        step = step.split("\n      - name:", 1)[0]
        script = step.split("        run: |\n", 1)[1]
        script = "\n".join(line[10:] for line in script.splitlines())
        info = dict(workflowName="Build all", event="push", status="completed",
                    conclusion="success", headSha="abc")
        info.update(overrides)
        env = dict(os.environ, RUN_ID="123", GITHUB_REPOSITORY="owner/repo",
                   TAG="v2.6.0-rme.10", WAIT_FOR_BUILD=str(wait).lower(),
                   MOCK_INFO=json.dumps(info), MOCK_SHA=tag_sha,
                   MOCK_REF=json.dumps({"object": {
                       "type": "tag" if annotated else "commit",
                       "sha": "annotated-object" if annotated else tag_sha}}))
        # Shell functions mock gh and sleep while exercising real bash/jq parsing.
        mocks = r'''
          gh() {
            if [ "$1 $2" = 'run view' ]; then
              if [ "$pending" = true ] && [ ! -f "$marker" ]; then
                printf '%s' "$MOCK_INFO" | jq '.status = "in_progress"'
              else
                printf '%s' "$MOCK_INFO"
              fi
            elif [[ "$2" == */git/ref/tags/* ]]; then
              printf '%s' "$MOCK_REF"
            elif [[ "$2" == */git/tags/* ]]; then
              printf '%s' "$MOCK_SHA"
            else
              echo 'Unexpected gh request' >&2; return 1
            fi
          }
          sleep() { touch "$marker"; }
        '''
        with tempfile.TemporaryDirectory() as tmp:
            env.update(marker=str(Path(tmp) / "waited"), pending=str(pending).lower())
            result = subprocess.run(["bash", "-c", mocks + script], env=env,
                                    capture_output=True, text=True, timeout=5)
        self.assertEqual(result.returncode == 0, accepted, result.stdout + result.stderr)

    def test_successful_matching_push(self):
        self.check_gate(accepted=True)

    def test_successful_manual_annotated_tag(self):
        self.check_gate(accepted=True, event="workflow_dispatch", annotated=True)

    def test_wrong_workflow(self):
        self.check_gate(workflowName="Check profiles")

    def test_pull_request(self):
        self.check_gate(event="pull_request")

    def test_wrong_commit(self):
        self.check_gate(tag_sha="old-commit")

    def test_failed_or_cancelled_build(self):
        for conclusion in ("failure", "cancelled", "skipped", ""):
            with self.subTest(conclusion=conclusion):
                self.check_gate(conclusion=conclusion)

    def test_attach_only_rejects_unfinished_build(self):
        self.check_gate(status="in_progress")

    def test_stable_waits_for_success(self):
        self.check_gate(accepted=True, wait=True, pending=True)

    def test_stable_waits_but_rejects_failure(self):
        self.check_gate(wait=True, pending=True, conclusion="failure")


if __name__ == "__main__":
    unittest.main()
