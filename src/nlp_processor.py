import json
from openai import OpenAI
from src.config import Config

class NLPProcessor:
    def __init__(self):
        self.client = OpenAI(base_url=Config.LM_STUDIO_API_BASE, api_key=Config.LM_STUDIO_API_KEY)

    def screen_paper(self, title, abstract, criteria_text):
        """
        Asks the LLM if the paper should be included based on the criteria.
        Returns: (included: bool, reason: str)
        """
        if not abstract or len(abstract) < 20:
             # Without abstract, hard to screen. 
             # Strategy: Include for manual review? Or exclude?
             # Let's mark as Undecided/Include to be safe.
             return True, "No abstract available - Included for manual check"

        prompt = f"""
        You are a research assistant screening papers for a systematic review.
        
        CRITERIA:
        {criteria_text}
        
        PAPER:
        Title: {title}
        Abstract: {abstract}
        
        TASK:
        Decide if this paper should be INCLUDED or EXCLUDED based on the criteria.
        Explain your reasoning briefly.
        
        OUTPUT FORMAT:
        Return a single JSON object:
        {{
          "included": true,
          "reason": "Brief explanation..."
        }}
        """
        
        try:
            response = self.client.chat.completions.create(
                model="local-model",
                messages=[
                    {"role": "system", "content": "You are a helpful research screening assistant. Output JSON only."},
                    {"role": "user", "content": prompt}
                ],
                temperature=0.0, # Low temp for deterministic logic
            )
            content = response.choices[0].message.content
            
            # Simple parsing
            # Clean md blocks
            content = content.replace("```json", "").replace("```", "").strip()
            
            # Attempt generic JSON search if strict parsing fails
            import re
            match = re.search(r'\{.*\}', content, re.DOTALL)
            if match:
                data = json.loads(match.group(0))
                return data.get("included", False), data.get("reason", "No reason provided")
            else:
                 return False, f"LLM Error (Format): {content[:50]}..."

        except Exception as e:
            return False, f"LLM Connection Error: {str(e)}"
