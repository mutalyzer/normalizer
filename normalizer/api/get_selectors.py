from flask_restx import Namespace, Resource, api, fields, reqparse

from normalizer.reference import get_reference_model, get_selectors_ids

ns = Namespace("/")


@ns.route("/get_selectors/<string:reference_id>")
class GetSelectors(Resource):
    def get(self, reference_id):
        """Retrieve available selectors for the provided reference."""
        reference_model = get_reference_model(reference_id)
        if reference_model:
            selectors = get_selectors_ids(reference_model["annotations"])
            return {"reference": reference_id, "selectors": selectors}
        return {"errors": [{"code": "ERETR"}]}
