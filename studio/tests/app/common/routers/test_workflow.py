workspace_id = "default"
unique_id = "0123"


def test_import(client):
    response = client.get(f"/workflow/reproduce/{workspace_id}/{unique_id}")
    data = response.json()

    assert response.status_code == 200
    assert isinstance(data, dict)
