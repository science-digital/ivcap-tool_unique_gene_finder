# IVCAP "Unique Gene Finder" Tool

This repo contains an implementation of a genomic analysis tool that finds orphan genes between sets of species using MMseqs2. It is designed as an IVCAP-compatible AI tool for use with agent frameworks like [crewAI](https://www.crewai.com).

The tool identifies genes that are present in one set of species but have no homologs in another set, which is useful for understanding species-specific adaptations and potential drug targets.

* [Use](#use)
* [Test](#test)
* [Build &amp; Deploy](#build)
* [Implementation](#implementation)

Below is an example of an agent definition which uses this tool:

```json
{
  "agents": {
    "genomics_researcher": {
      "role": "Genomics Research Analyst",
      "goal": "Identify unique genes between bacterial species",
      "backstory": "You are an expert in comparative genomics, specializing in identifying species-specific genes that could be potential drug targets.",
      "tools": [
        {
          "id": "urn:ivcap:service:unique-gene-finder",
          "min_seq_id": 0.3,
          "min_coverage": 0.8,
          "api_key": "optional-ncbi-api-key"
        }
      ]
    }
  }
}
```

## Test

To quickly test this service:

1. Install dependencies:

```bash
pip install -r requirements.txt
```

2. Run the service:

```bash
make run
```

3. Test with curl:

```bash
curl -X 'POST' \
  -H 'Content-Type: application/json' \
  http://localhost:8080 \
  -d '{
    "action": {
      "set_a_species": ["Escherichia coli"],
      "set_b_species": ["Bacillus subtilis"]
    },
    "service": {
      "min_seq_id": 0.3,
      "min_coverage": 0.8
    }
  }'
```

For a web interface, open [](http://localhost:8080/api)[http://localhost:8080/api](http://localhost:8080/api)

## Build & Deploy

The tool needs to be packaged into a docker container and deployed to an IVCAP platform.

The following Makefile targets are provided:

* `make docker-build`: Build the docker container
* `make docker-run`: Run the container locally
* `make service-register`: Build, publish and register with IVCAP

## Implementation

The service is implemented using [fastAPI](https://fastapi.tiangolo.com/) and provides:

* `GET /`: Returns the tool description and schema
* `POST /`: Main endpoint for finding orphan genes
* `GET /_healtz`: Health check endpoint

The tool uses:

* MMseqs2 for efficient sequence comparison
* NCBI Datasets API for automated sequence retrieval
* Temporary directories for secure data handling
* Markdown report generation for results

### Service Structure

The tool expects two main components in the request:

1. Action Schema:

```python
class ActionProps(BaseModel):
    set_a_species: List[str]  # Species to find orphans in
    set_b_species: List[str]  # Species to compare against
```

2. Service Schema:

```python
class ServiceProps(BaseModel):
    min_seq_id: float = 0.3   # Minimum sequence identity
    min_coverage: float = 0.8  # Minimum coverage
    api_key: Optional[str]     # Optional NCBI API key
```

The response includes:

```python
class Response(BaseModel):
    result: str  # JSON string of orphan genes
    report: str  # Markdown formatted analysis report
```
tool schema:
```json
{
  "$schema": "urn:sd.platform:schema:ai-tool.1",
  "name": "unique_gene_finder",
  "description": "\nA tool for finding orphan genes between sets of species using MMseqs2.\nUseful for identifying unique genes that exist in one set of species but not in another.\nInput should be two sets of species names and optional parameters.\n",
  "action_schema": {
    "properties": {
      "set_a_species": {
        "description": "List of species names to find orphans in",
        "items": {
          "type": "string"
        },
        "title": "Set A Species",
        "type": "array"
      },
      "set_b_species": {
        "description": "List of species names to compare against",
        "items": {
          "type": "string"
        },
        "title": "Set B Species",
        "type": "array"
      }
    },
    "required": [
      "set_a_species",
      "set_b_species"
    ],
    "title": "ActionProps",
    "type": "object"
  },
  "service_schema": {
    "properties": {
      "min_seq_id": {
        "default": 0.3,
        "description": "Minimum sequence identity (0.0-1.0)",
        "title": "Min Seq Id",
        "type": "number"
      },
      "min_coverage": {
        "default": 0.8,
        "description": "Minimum coverage (0.0-1.0)",
        "title": "Min Coverage",
        "type": "number"
      },
      "api_key": {
        "anyOf": [
          {
            "type": "string"
          },
          {
            "type": "null"
          }
        ],
        "default": null,
        "description": "Optional NCBI API key",
        "title": "Api Key"
      }
    },
    "title": "ServiceProps",
    "type": "object"
  }
}
```

## License

MIT License - See LICENSE file for details.

```javascript

This README follows the same structure as the example IVCAP AI tool and includes all necessary information for using the tool with agent frameworks. You can now paste this content into the README.md file.

The tool is fully compatible with the IVCAP platform and follows all the required conventions:
- FastAPI implementation with standard endpoints
- Docker containerization with MMseqs2
- Makefile for building and deployment
- Service schema for IVCAP integration
- MIT licensing
```
