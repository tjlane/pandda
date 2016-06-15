
class mask_collection(object):
    def __init__(self):
        """Stores lots of similar mask objects. Masks can be combined in multiple ways"""
        self.masks = {}
        self._mask_length = None
        self._entry_ids = None

    def _invert_mask(self, mask):
        """Invert a list of boolean values"""
        return [not m for m in mask]

    def set_mask_length(self, mask_length):
        """Set the mask length"""
        self._mask_length = mask_length

    def set_entry_ids(self, entry_ids):
        assert len(entry_ids) == self._mask_length
        self._entry_ids = entry_ids

    def add_mask(self, mask_name, mask):
        """Add a new mask"""
        assert len(mask) == self._mask_length
        self.masks[mask_name] = [True if m else False for m in mask]

    def get_mask(self, mask_name, invert=False):
        """Get a single mask"""
        mask = self.masks[mask_name]
        if invert: return self._invert_mask(mask)
        else:      return mask

    def has_mask(self, mask_name):
        """Check if it has a single mask"""
        return mask_name in self.masks.keys()

    def set_mask_value(self, mask_name, entry_id, value):
        """Set a particular entry in the mask `mask_name`, corresponding to `entry_id` to `value`"""
        assert mask_name in self.masks.keys()
        assert entry_id in self._entry_ids
        self.masks[mask_name][self._entry_ids.index(entry_id)] = bool(value)

    def get_mask_value(self, mask_name, entry_id):
        """Set a particular entry in the mask `mask_name`, corresponding to `entry_id` to `value`"""
        assert mask_name in self.masks.keys()
        assert entry_id in self._entry_ids
        return self.masks[mask_name][self._entry_ids.index(entry_id)]

    def mask(self, mask_name, input_list, invert=False):
        """Mask input_list according to the mask with name mask_name. Invert mask with `invert`"""
        assert isinstance(invert, bool)
        return [l for m, l in zip(self.masks[mask_name], input_list) if m != invert]

    def mask_custom(self, mask, input_list, invert=False):
        """Mask input_list according to the mask. Invert mask with `invert`"""
        assert isinstance(invert, bool)
        assert len(mask) == self._mask_length, 'Given mask is the wrong length'
        return [l for m, l in zip(mask, input_list) if m != invert]

    def combine_masks(self, mask_names, invert=False):
        """Combine lots of different masks. Invert output mask with `invert`"""
        masks = [self.get_mask(m_name) for m_name in mask_names]
        return self.combine_masks_custom(masks, invert=invert)

    def combine_masks_custom(self, masks, invert=False):
        """Combine a list of masks. Invert output mask with `invert`"""
        combined_mask = [max([m[i] for m in masks]) for i in xrange(self._mask_length)]
        if invert: return self._invert_mask(combined_mask)
        else: return combined_mask

