"""Thin client for the Open Targets Platform GraphQL API (public, no auth/license gate).

Used for target validation: genetic evidence, safety liabilities, paralog
family (selectivity risk), baseline expression/localization, tractability.
"""
import requests

API_URL = "https://api.platform.opentargets.org/api/v4/graphql"

SEARCH_QUERY = """
query Search($q: String!) {
  search(queryString: $q, entityNames: ["target"]) {
    hits { id name entity }
  }
}
"""

TARGET_PROFILE_QUERY = """
query TargetProfile($id: String!) {
  target(ensemblId: $id) {
    id
    approvedSymbol
    approvedName
    biotype
    functionDescriptions
    targetClass { id label }
    tractability { label modality value }
    geneticConstraint { constraintType exp obs score oe oeLower oeUpper }
    isEssential
    safetyLiabilities {
      event eventId datasource literature url
      effects { dosing direction }
      biosamples { tissueLabel cellLabel }
      studies { name type description }
    }
    homologues {
      speciesName targetGeneSymbol targetGeneId homologyType
      queryPercentageIdentity targetPercentageIdentity isHighConfidence
    }
    subcellularLocations { location termSL source }
    baselineExpression {
      count
      rows { tissueBiosample { biosampleName } median unit }
    }
    associatedDiseases {
      count
      rows { disease { id name } score datatypeScores { id score } }
    }
  }
}
"""


def resolve_ensembl_id(gene_symbol: str) -> str:
    resp = requests.post(
        API_URL, json={"query": SEARCH_QUERY, "variables": {"q": gene_symbol}}, timeout=30
    )
    resp.raise_for_status()
    hits = resp.json()["data"]["search"]["hits"]
    for hit in hits:
        if hit["name"].upper() == gene_symbol.upper():
            return hit["id"]
    raise ValueError(f"No exact target match for {gene_symbol!r} (got: {[h['name'] for h in hits]})")


def fetch_target_profile(ensembl_id: str) -> dict:
    resp = requests.post(
        API_URL, json={"query": TARGET_PROFILE_QUERY, "variables": {"id": ensembl_id}}, timeout=30
    )
    resp.raise_for_status()
    payload = resp.json()
    if "errors" in payload:
        raise RuntimeError(payload["errors"])
    return payload["data"]["target"]
